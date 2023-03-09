#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/EmptyGroupDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/CRC16.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"

#include "DataFormats/HGCalDigi/interface/HGCalRawDataEmulatorInfo.h"
#include "EventFilter/HGCalRawToDigi/interface/HGCalFrameGenerator.h"
#include "EventFilter/HGCalRawToDigi/interface/TBTreeReader.h"

class HGCalSlinkEmulator : public edm::stream::EDProducer<> {
public:
  explicit HGCalSlinkEmulator(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const unsigned int fed_id_;

  const bool store_emul_info_;
  const bool store_fed_header_trailer_;

  std::unique_ptr<hgcal::econd::Emulator> emulator_;

  edm::Service<edm::RandomNumberGenerator> rng_;
  edm::EDPutTokenT<FEDRawDataCollection> fedRawToken_;
  edm::EDPutTokenT<HGCalSlinkEmulatorInfo> fedEmulInfoToken_;
  hgcal::HGCalFrameGenerator frame_gen_;
};

HGCalSlinkEmulator::HGCalSlinkEmulator(const edm::ParameterSet& iConfig)
    : fed_id_(iConfig.getParameter<unsigned int>("fedId")),
      store_emul_info_(iConfig.getParameter<bool>("storeEmulatorInfo")),
      store_fed_header_trailer_(iConfig.getParameter<bool>("fedHeaderTrailer")),
      frame_gen_(iConfig) {
  // figure out which emulator is to be used
  const auto& emul_type = iConfig.getParameter<std::string>("emulatorType");
  if (emul_type == "empty")
    emulator_ = std::make_unique<hgcal::econd::EmptyEmulator>(frame_gen_.econdParams().num_channels);
  else if (emul_type == "trivial")
    emulator_ = std::make_unique<hgcal::econd::TrivialEmulator>(
        frame_gen_.econdParams().num_channels, iConfig.getUntrackedParameter<std::vector<unsigned int>>("channels"));
  else if (emul_type == "tbTree")
    emulator_ =
        std::make_unique<hgcal::econd::TBTreeReader>(iConfig.getUntrackedParameter<std::string>("treeName"),
                                                     iConfig.getUntrackedParameter<std::vector<std::string>>("inputs"),
                                                     frame_gen_.econdParams().num_channels);
  else
    throw cms::Exception("HGCalSlinkEmulator") << "Invalid emulator type chosen: '" << emul_type << "'.";

  // ensure the random number generator service is present in configuration
  if (!rng_.isAvailable())
    throw cms::Exception("HGCalSlinkEmulator") << "The HGCalSlinkEmulator module requires the "
                                                  "RandomNumberGeneratorService,\n"
                                                  "which appears to be absent. Please add that service to your "
                                                  "configuration\n"
                                                  "or remove the modules that require it.";

  fedRawToken_ = produces<FEDRawDataCollection>();
  if (store_emul_info_)
    fedEmulInfoToken_ = produces<HGCalSlinkEmulatorInfo>();
}

void HGCalSlinkEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const auto reader_evt = emulator_->next();
  const auto event_id = std::get<0>(reader_evt.first), bx_id = std::get<1>(reader_evt.first);

  frame_gen_.setRandomEngine(rng_->getEngine(iEvent.streamID()));

  // build the S-link payload
  auto slink_event = frame_gen_.produceSlinkEvent(fed_id_, reader_evt);
  const auto slink_event_size = slink_event.size() * sizeof(slink_event.at(0));

  // compute the total S-link payload size
  size_t total_event_size = slink_event_size;
  if (store_fed_header_trailer_)
    total_event_size += FEDHeader::length + FEDTrailer::length;

  // fill the output FED raw data collection
  FEDRawDataCollection raw_data;
  auto& fed_data = raw_data.FEDData(fed_id_);
  fed_data.resize(total_event_size);
  auto* ptr = fed_data.data();

  if (store_fed_header_trailer_) {
    int trg_type = 0;  //FIXME
    // compose 2*32-bit FED header word
    FEDHeader::set(ptr, trg_type, event_id, bx_id, fed_id_);
    LogDebug("HGCalSlinkEmulator").log([&](auto& log) {
      const FEDHeader hdr(ptr);
      log << "FED header: lvl1ID=" << hdr.lvl1ID() << ", bxID=" << hdr.bxID() << ", source ID=" << hdr.sourceID()
          << ".";
    });
    ptr += FEDHeader::length;
  }

  // insert ECON-D payload
  LogDebug("HGCalSlinkEmulator") << "Will write " << slink_event.size() << " 64-bit words = " << slink_event_size
                                 << " 8-bit words.";
  /*for (const auto& word : slink_event)
    for (ssize_t i = 7; i >= 0; --i) {
      *ptr = (word >> i * 8) & 0xff;
      ++ptr;
    }*/
  std::memcpy(ptr, slink_event.data(), slink_event_size);
  ptr += slink_event_size;

  if (store_fed_header_trailer_) {
    // compose 2*32-bit FED trailer word
    FEDTrailer::set(ptr,
                    slink_event.size() + 2,
                    evf::compute_crc(reinterpret_cast<uint8_t*>(slink_event.data()), slink_event_size),
                    0,
                    0);
    LogDebug("HGCalSlinkEmulator").log([&](auto& log) {
      const FEDTrailer trl(ptr);
      log << "FED trailer: fragment length: " << trl.fragmentLength() << ", CRC=0x" << std::hex << trl.crc() << std::dec
          << ", status: " << trl.evtStatus() << ".";
    });
    ptr += FEDTrailer::length;
  }

  iEvent.emplace(fedRawToken_, std::move(raw_data));

  // store the emulation information if requested
  if (store_emul_info_) {
    auto emul_info = frame_gen_.lastSlinkEmulatedInfo();
    iEvent.emplace(fedEmulInfoToken_, std::move(emul_info));
  }
}

void HGCalSlinkEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  auto desc = hgcal::HGCalFrameGenerator::description();
  desc.ifValue(edm::ParameterDescription<std::string>("emulatorType", "empty", true),
               // empty payload emulator
               "empty" >> edm::EmptyGroupDescription() or
                   // trivial emulator
                   "trivial" >> edm::ParameterDescription<std::vector<unsigned int>>("channels", {}, false) or
                   // test beam tree content
                   "tbTree" >> (edm::ParameterDescription<std::string>("treeName", "unpacker_data/hgcroc", false) and
                                edm::ParameterDescription<std::vector<std::string>>("inputs", {}, false)));
  desc.add<unsigned int>("fedId", 0)->setComment("FED number delivering the emulated frames");
  desc.add<bool>("fedHeaderTrailer", false)->setComment("also add FED header/trailer info");
  desc.add<bool>("storeEmulatorInfo", true)
      ->setComment("also append a 'truth' auxiliary info to the output event content");
  descriptions.add("hgcalEmulatedSlinkRawData", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalSlinkEmulator);
