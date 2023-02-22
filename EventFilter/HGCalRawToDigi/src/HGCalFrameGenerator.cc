#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/CRC16.h"

#include "EventFilter/HGCalRawToDigi/interface/HGCalFrameGenerator.h"
#include "EventFilter/HGCalRawToDigi/interface/HGCalRawDataPackingTools.h"

#include "CLHEP/Random/RandFlat.h"

namespace hgcal {
  HGCalFrameGenerator::HGCalFrameGenerator(const edm::ParameterSet& iConfig)
      : pass_through_(iConfig.getParameter<bool>("passThroughMode")),
        expected_mode_(iConfig.getParameter<bool>("expectedMode")),
        matching_ebo_numbers_(iConfig.getParameter<bool>("matchingEBOnumbers")),
        bo_truncated_(iConfig.getParameter<bool>("bufferOverflowTruncated")) {
    const auto econd_params = iConfig.getParameter<edm::ParameterSet>("econdParams");
    econd_.chan_surv_prob = econd_params.getParameter<double>("channelSurv");
    econd_.enabled_channels = econd_params.getParameter<std::vector<unsigned int> >("enabledChannels");
    econd_.header_marker = econd_params.getParameter<unsigned int>("headerMarker");
    econd_.num_channels = econd_params.getParameter<unsigned int>("numChannels");
    econd_.bitO_error_prob = econd_params.getParameter<double>("bitOError");
    econd_.bitB_error_prob = econd_params.getParameter<double>("bitBError");
    econd_.bitE_error_prob = econd_params.getParameter<double>("bitEError");
    econd_.bitT_error_prob = econd_params.getParameter<double>("bitTError");
    econd_.bitH_error_prob = econd_params.getParameter<double>("bitHError");
    econd_.bitS_error_prob = econd_params.getParameter<double>("bitSError");

    const auto slink_params = iConfig.getParameter<edm::ParameterSet>("slinkParams");
    slink_.num_econds = slink_params.getParameter<unsigned int>("numECONDs");
  }

  edm::ParameterSetDescription HGCalFrameGenerator::description() {
    edm::ParameterSetDescription desc;

    desc.add<bool>("passThroughMode", true)->setComment("ECON-D in pass-through mode?");
    desc.add<bool>("expectedMode", false)->setComment("is an Event HDR/TRL expected to be received from the HGCROCs?");
    desc.add<bool>("matchingEBOnumbers", false)
        ->setComment(
            "is the transmitted E/B/O (according to mode selected by user) matching the E/B/O value in the ECON-D L1A "
            "FIFO?");
    desc.add<bool>("bufferOverflowTruncated", false)->setComment("is the packet truncated for buffer overflow?");

    edm::ParameterSetDescription econd_desc;
    econd_desc.add<double>("channelSurv", 1.);
    econd_desc.add<std::vector<unsigned int> >("enabledChannels", {})
        ->setComment("list of channels to be enabled in readout");
    econd_desc.add<unsigned int>("headerMarker", 0x154)
        ->setComment("9b programmable pattern; default is '0xAA' + '0b0'");
    econd_desc.add<unsigned int>("numChannels", 37)->setComment("number of channels managed in ECON-D");
    econd_desc.add<double>("bitOError", 0.);
    econd_desc.add<double>("bitBError", 0.);
    econd_desc.add<double>("bitEError", 0.);
    econd_desc.add<double>("bitTError", 0.);
    econd_desc.add<double>("bitHError", 0.);
    econd_desc.add<double>("bitSError", 0.);
    desc.add<edm::ParameterSetDescription>("econdParams", econd_desc);

    edm::ParameterSetDescription slink_desc;
    slink_desc.add<unsigned int>("numECONDs", 7);
    desc.add<edm::ParameterSetDescription>("slinkParams", slink_desc);

    return desc;
  }

  void HGCalFrameGenerator::setRandomEngine(CLHEP::HepRandomEngine& rng) { rng_ = &rng; }

  HGCalFrameGenerator::HeaderBits HGCalFrameGenerator::generateStatusBits() const {
    HeaderBits header_bits;
    // first sample on header status bits
    header_bits.bitO = CLHEP::RandFlat::shoot(rng_) >= econd_.bitO_error_prob;
    header_bits.bitB = CLHEP::RandFlat::shoot(rng_) >= econd_.bitB_error_prob;
    header_bits.bitE = CLHEP::RandFlat::shoot(rng_) >= econd_.bitE_error_prob;
    header_bits.bitT = CLHEP::RandFlat::shoot(rng_) >= econd_.bitT_error_prob;
    header_bits.bitH = CLHEP::RandFlat::shoot(rng_) >= econd_.bitH_error_prob;
    header_bits.bitS = CLHEP::RandFlat::shoot(rng_) >= econd_.bitS_error_prob;
    return header_bits;
  }

  std::vector<bool> HGCalFrameGenerator::generateEnabledChannels(uint64_t& ch_en) const {
    std::vector<bool> chmap(econd_.num_channels, false);
    ch_en = 0ull;  // reset the list of channels enabled
    for (size_t i = 0; i < chmap.size(); i++) {
      // randomly choosing the channels to be shot at
      chmap[i] = (econd_.enabled_channels.empty() ||
                  (std::find(econd_.enabled_channels.begin(), econd_.enabled_channels.end(), i) !=
                   econd_.enabled_channels.end())) &&
                 CLHEP::RandFlat::shoot(rng_) <= econd_.chan_surv_prob;
      ch_en += (chmap[i] << i);
    }
    ch_en &= ((1 << (econd_.num_channels + 1)) - 1);  // mask only (econd_.num_channels) LSBs
    return chmap;
  }

  std::vector<uint32_t> HGCalFrameGenerator::generateERxData(const econd::ERxEvent& event,
                                                             std::vector<uint64_t>& enabled_channels) const {
    enabled_channels.clear();  // reset the list of channels enabled

    std::vector<uint32_t> erxData;
    uint64_t ch_en;  // list of channels enabled
    for (const auto& jt : event) {
      auto chmap = generateEnabledChannels(ch_en);
      enabled_channels.emplace_back(ch_en);

      auto erxHeader = econd::eRxSubPacketHeader(0, 0, false, jt.second.cm0, jt.second.cm1, chmap);
      erxData.insert(erxData.end(), erxHeader.begin(), erxHeader.end());
      if (jt.second.adc.size() < econd_.num_channels) {
        edm::LogVerbatim("HGCalFrameGenerator:generateERxData")
            << "Data multiplicity to low (" << jt.second.adc.size() << ") to emulate " << econd_.num_channels
            << " ECON-D channel(s).";
        continue;
      }
      for (size_t i = 0; i < econd_.num_channels; i++) {
        if (!chmap.at(i))
          continue;
        uint8_t msb = 32;
        auto chData = econd::addChannelData(msb,
                                            jt.second.tctp.at(i),
                                            jt.second.adc.at(i),
                                            jt.second.tot.at(i),
                                            jt.second.adcm.at(i),
                                            jt.second.toa.at(i),
                                            true,
                                            true,
                                            true,
                                            true);
        erxData.insert(erxData.end(), chData.begin(), chData.end());
      }
    }
    return erxData;
  }

  std::vector<uint32_t> HGCalFrameGenerator::produceECONEvent(const econd::ECONDEvent& event) const {
    std::vector<uint64_t> enabled_channels;
    auto header_bits = generateStatusBits();
    auto econd_event = generateERxData(event.second, enabled_channels);
    LogDebug("HGCalFrameGenerator") << econd_event.size() << " word(s) of eRx payloads inserted.";

    last_econd_emul_info_.clear();
    // as ECON-D event content was just created, only prepend packet header at this stage
    const auto econd_header = hgcal::econd::eventPacketHeader(
        econd_.header_marker,
        econd_event.size() + 1,
        pass_through_,
        expected_mode_,
        // HGCROC Event reco status across all active eRxE-B-O:
        // FIXME check endianness of these two numbers
        (header_bits.bitH << 1) | header_bits.bitT,                            // HDR/TRL numbers
        (header_bits.bitE << 2) | (header_bits.bitB << 1) | header_bits.bitO,  // Event/BX/Orbit numbers
        matching_ebo_numbers_,
        bo_truncated_,
        0,                         // Hamming for event header
        std::get<1>(event.first),  // BX
        std::get<0>(event.first),  // event id (L1A)
        std::get<2>(event.first),  // orbit
        header_bits.bitS,          // OR of "Stat" bits for all active eRx
        0,
        0  // CRC
    );
    econd_event.insert(econd_event.begin(), econd_header.begin(), econd_header.end());
    LogDebug("HGCalFrameGenerator") << econd_header.size()
                                    << " word(s) of event packet header prepend. New size of ECON frame: "
                                    << econd_event.size();

    //econd_event.push_back(computeCRC(econd_header));

    return econd_event;
  }

  static std::vector<uint64_t> to64bit(const std::vector<uint32_t>& in) {
    std::vector<uint64_t> out;
    for (size_t i = 0; i < in.size(); i += 2) {
      uint64_t word1 = in.at(i), word2 = (i + 1 < in.size()) ? in.at(i + 1) : 0ul;
      out.emplace_back(((word2 & 0xffffffff) << 32) | (word1 & 0xffffffff));
    }
    return out;
  }

  std::vector<uint64_t> HGCalFrameGenerator::produceSlinkEvent(uint32_t fed_id,
                                                               const econd::ECONDEvent& econd_event) const {
    std::vector<uint64_t> slink_event;

    const auto& eid = econd_event.first;
    uint64_t event_id = std::get<0>(eid), bx_id = std::get<1>(eid), orbit_id = std::get<2>(eid);

    // build the S-link header words
    const uint32_t content_id = 0;
    const uint8_t r6 = 0, r8 = 0, v = 0, boe = 0;

    const auto slink_header = to64bit(backend::buildSlinkHeader(boe, v, r8, event_id, r6, content_id, fed_id));

    uint64_t l1a_header{0ull};
    l1a_header |= ((orbit_id & 0xf) << 36);
    l1a_header |= ((event_id & 0x3f) << 40);
    l1a_header |= ((bx_id & 0xfff) << 46);

    last_slink_emul_info_.clear();
    for (size_t i = 0; i < max_num_econds_; ++i) {
      if (i < slink_.num_econds) {  // active ECON-D
        auto econd_evt = to64bit(produceECONEvent(econd_event));
        const auto& econd_evt_info = lastECONDEmulatedInfo();
        slink_event.insert(slink_event.end(), econd_evt.begin(), econd_evt.end());
        uint8_t econd_packet_status = ECONDPacketStatus::Normal;  //FIXME
        l1a_header |= ((econd_packet_status & 0x7) << (3 * i));
        last_slink_emul_info_.addECONDEmulatedInfo(econd_evt_info);
      } else {  // inactive ECON-D
        l1a_header |= (ECONDPacketStatus::InactiveECOND << (3 * i));
      }
    }
    uint64_t l1a_header_reversed = ((l1a_header & 0xffffffff) << 32) | ((l1a_header >> 32) & 0xffffffff);
    slink_event.insert(slink_event.begin(), l1a_header_reversed);                       // prepend L1A header
    slink_event.insert(slink_event.begin(), slink_header.begin(), slink_header.end());  // prepend S-link header

    // build the S-link trailer words
    uint8_t eoe = 0, daqcrc = 0, trailer_r = 0;
    uint64_t event_length = slink_event.size() - slink_header.size() - 1;
    uint32_t crc = 0, status = 0;
    const auto slink_trailer =
        to64bit(backend::buildSlinkTrailer(eoe, daqcrc, trailer_r, event_length, bx_id, orbit_id, crc, status));
    slink_event.insert(slink_event.end(), slink_trailer.begin(), slink_trailer.end());

    return slink_event;
  }

  uint8_t HGCalFrameGenerator::computeCRC(const std::vector<uint32_t>& event_header) const {
    uint8_t crc = 0;  //FIXME 8-bit Bluetooth CRC
    return crc;
  }
}  // namespace hgcal
