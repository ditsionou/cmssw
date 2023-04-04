/****************************************************************************
 *
 * This is a part of HGCAL offline software.
 * Authors:
 *   Laurent Forthomme
 *
 ****************************************************************************/

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/CRC16.h"

#include "EventFilter/HGCalRawToDigi/interface/HGCalECONDEmulator.h"
#include "EventFilter/HGCalRawToDigi/interface/HGCalFrameGenerator.h"
#include "EventFilter/HGCalRawToDigi/interface/HGCalRawDataPackingTools.h"

#include "CLHEP/Random/RandFlat.h"
#include <iomanip>

namespace hgcal {
  // utilities
  template <typename T>
  void printWords(edm::MessageSender& os, const std::string& name, const std::vector<T> vec) {
    os << "Dump of the '" << name << "' words:\n";
    for (size_t i = 0; i < vec.size(); ++i)
      os << std::dec << std::setfill(' ') << std::setw(4) << i << ": 0x" << std::hex << std::setfill('0')
         << std::setw(sizeof(T) * 2) << vec.at(i) << "\n";
  }

  static std::vector<uint64_t> to64bit(const std::vector<uint32_t>& in) {
    std::vector<uint64_t> out;
    for (size_t i = 0; i < in.size(); i += 2) {
      uint64_t word1 = (i < in.size()) ? in.at(i) : 0ul, word2 = (i + 1 < in.size()) ? in.at(i + 1) : 0ul;
      out.emplace_back(((word2 & 0xffffffff) << 32) | (word1 & 0xffffffff));
    }
    return out;
  }

  static std::vector<uint64_t> to128bit(const std::vector<uint32_t>& in) {
    std::vector<uint64_t> out;
    for (size_t i = 0; i < in.size(); i += 4) {
      uint64_t word1 = (i < in.size()) ? in.at(i) : 0ul, word2 = (i + 1 < in.size()) ? in.at(i + 1) : 0ul,
               word3 = (i + 2 < in.size()) ? in.at(i + 2) : 0ul, word4 = (i + 3 < in.size()) ? in.at(i + 3) : 0ul;
      out.emplace_back(((word2 & 0xffffffff) << 32) | (word1 & 0xffffffff));
      out.emplace_back(((word4 & 0xffffffff) << 32) | (word3 & 0xffffffff));
    }
    return out;
  }

  HGCalFrameGenerator::HGCalFrameGenerator(const edm::ParameterSet& iConfig) {
    const auto slink_params = iConfig.getParameter<edm::ParameterSet>("slinkParams");

    size_t econd_id = 0;
    std::vector<unsigned int> active_econds;
    for (const auto& econd : slink_params.getParameter<std::vector<edm::ParameterSet> >("ECONDs")) {
      // a bit of user input validation
      if (active_econds.size() > max_num_econds_)
        throw cms::Exception("HGCalFrameGenerator")
            << "Too many active ECON-D set: " << slink_.active_econds.size() << " > " << max_num_econds_ << ".";
      if (econd_id >= max_num_econds_)
        throw cms::Exception("HGCalFrameGenerator")
            << "Invalid ECON-D identifier: " << econd_id << " >= " << max_num_econds_ << ".";
      if (econd.getParameter<bool>("active"))
        active_econds.emplace_back(econd_id);

      econd_params_.insert(std::make_pair(econd_id, econd::EmulatorParameters(econd)));
      ++econd_id;
    }

    slink_ = SlinkParameters{.active_econds = active_econds,
                             .boe_marker = slink_params.getParameter<unsigned int>("boeMarker"),
                             .eoe_marker = slink_params.getParameter<unsigned int>("eoeMarker"),
                             .format_version = slink_params.getParameter<unsigned int>("formatVersion"),
                             .num_capture_blocks = slink_params.getParameter<unsigned int>("numCaptureBlocks"),
                             .store_header_trailer = slink_params.getParameter<bool>("storeHeaderTrailer")};
  }

  const econd::EmulatorParameters& HGCalFrameGenerator::econdParams(unsigned int econd_id) const {
    if (econd_params_.count(econd_id) == 0)
      throw cms::Exception("HGCalFrameGenerator") << "ECON-D with id=" << econd_id << " was not configured.";
    return econd_params_.at(econd_id);
  }

  edm::ParameterSetDescription HGCalFrameGenerator::description() {
    edm::ParameterSetDescription desc;

    std::vector<edm::ParameterSet> econds_psets;
    for (size_t i = 0; i < 7; ++i)
      econds_psets.emplace_back();

    edm::ParameterSetDescription slink_desc;
    slink_desc.addVPSet("ECONDs", econd::EmulatorParameters::description(), econds_psets)
        ->setComment("list of active ECON-Ds in S-link");
    slink_desc.add<unsigned int>("boeMarker", 0x55);
    slink_desc.add<unsigned int>("eoeMarker", 0xaa);
    slink_desc.add<unsigned int>("formatVersion", 3);
    slink_desc.add<unsigned int>("numCaptureBlocks", 1)
        ->setComment("number of capture blocks to emulate per S-link payload");
    slink_desc.add<bool>("storeHeaderTrailer", true)->setComment("also store the S-link header and trailer words");
    desc.add<edm::ParameterSetDescription>("slinkParams", slink_desc);

    return desc;
  }

  void HGCalFrameGenerator::setRandomEngine(CLHEP::HepRandomEngine& rng) { rng_ = &rng; }

  void HGCalFrameGenerator::setEmulator(econd::Emulator& emul) { emul_ = &emul; }

  HGCalFrameGenerator::HeaderBits HGCalFrameGenerator::generateStatusBits(unsigned int econd_id) const {
    if (!rng_)
      throw cms::Exception("HGCalFrameGenerator") << "Random number generator not initialised.";
    const auto& econd_params = econdParams(econd_id);
    HeaderBits header_bits;
    // first sample on header status bits
    header_bits.bitO = CLHEP::RandFlat::shoot(rng_) < econd_params.error_prob.bitO;
    header_bits.bitB = CLHEP::RandFlat::shoot(rng_) < econd_params.error_prob.bitB;
    header_bits.bitE = CLHEP::RandFlat::shoot(rng_) < econd_params.error_prob.bitE;
    header_bits.bitT = CLHEP::RandFlat::shoot(rng_) < econd_params.error_prob.bitT;
    header_bits.bitH = CLHEP::RandFlat::shoot(rng_) < econd_params.error_prob.bitH;
    header_bits.bitS = CLHEP::RandFlat::shoot(rng_) < econd_params.error_prob.bitS;
    return header_bits;
  }

  std::vector<bool> HGCalFrameGenerator::generateEnabledChannels(unsigned int econd_id) const {
    const auto& econd_params = econdParams(econd_id);
    std::vector<bool> chmap(econd_params.num_channels_per_erx, false);
    for (size_t i = 0; i < chmap.size(); i++)
      // randomly choosing the channels to be shot at
      chmap[i] = CLHEP::RandFlat::shoot(rng_) <= econd_params.chan_surv_prob;
    return chmap;
  }

  std::vector<uint32_t> HGCalFrameGenerator::generateERxData(unsigned int econd_id,
                                                             const econd::ERxInput& event) const {
    const auto& econd_params = econdParams(econd_id);
    std::vector<uint32_t> erx_data;
    for (const auto& jt : event) {  // one per eRx
      auto chmap =
          generateEnabledChannels(econd_id);  // generate a list of probable channels to be filled with emulated content

      // insert eRx header (common mode, channels map, ...)
      uint8_t stat = 0b111 /*possibly emulate*/, hamming_check = 0;
      bool bit_e = false;
      auto erx_header = econd::eRxSubPacketHeader(stat, hamming_check, bit_e, jt.second.cm0, jt.second.cm1, chmap);
      erx_data.insert(erx_data.end(), erx_header.begin(), erx_header.end());
      if (jt.second.adc.size() < econd_params.num_channels_per_erx) {
        edm::LogWarning("HGCalFrameGenerator:generateERxData")
            << "Data multiplicity to low (" << jt.second.adc.size() << ") to emulate "
            << econd_params.num_channels_per_erx << " ECON-D channel(s).";
        continue;
      }
      // insert eRx payloads (integrating all readout channels)
      const auto erx_chan_data = econd::produceERxData(chmap,
                                                       jt.second,
                                                       true,  // passZS
                                                       true,  // passZSm1
                                                       true,  // hasToA
                                                       econd_params.characterisation_mode);
      erx_data.insert(erx_data.end(), erx_chan_data.begin(), erx_chan_data.end());
    }
    LogDebug("HGCalFrameGenerator").log([&erx_data](auto& log) { printWords(log, "erx", erx_data); });
    return erx_data;
  }

  std::vector<uint32_t> HGCalFrameGenerator::produceECONEvent(unsigned int econd_id) const {
    if (!emul_)
      throw cms::Exception("HGCalFrameGenerator")
          << "ECON-D emulator was not properly set to the frame generator. Please ensure you are calling the "
             "HGCalFrameGenerator::setEmulator method.";

    const auto event = emul_->next();
    last_emul_event_ = event;
    const auto& econd_params = econdParams(econd_id);
    auto header_bits = generateStatusBits(econd_id);
    auto econd_event = generateERxData(econd_id, event.second);
    LogDebug("HGCalFrameGenerator") << econd_event.size() << " word(s) of eRx payloads inserted.";

    last_econd_emul_info_.clear();

    // as ECON-D event content was just created, only prepend packet header at this stage
    auto econd_header =
        econd::eventPacketHeader(econd_params.header_marker,
                                 econd_event.size() + 1 /*CRC*/,
                                 econd_params.passthrough_mode,
                                 econd_params.expected_mode,
                                 // HGCROC Event reco status across all active eRx E-B-O:
                                 // (check endianness of these two numbers)
                                 (header_bits.bitH & 0x1) << 1 | (header_bits.bitT & 0x1),  // HDR/TRL numbers
                                 (header_bits.bitE & 0x1) << 2 | (header_bits.bitB & 0x1) << 1 |
                                     (header_bits.bitO & 0x1),  // Event/BX/Orbit numbers
                                 econd_params.matching_ebo_numbers,
                                 econd_params.bo_truncated,
                                 0,                         // Hamming for event header
                                 std::get<1>(event.first),  // BX
                                 std::get<0>(event.first),  // event id (L1A)
                                 std::get<2>(event.first),  // orbit
                                 header_bits.bitS,          // OR of "Stat" bits for all active eRx
                                 0);
    LogDebug("HGCalFrameGenerator").log([&econd_header](auto& log) { printWords(log, "econ-d header", econd_header); });
    econd_event.insert(econd_event.begin(), econd_header.begin(), econd_header.end());
    LogDebug("HGCalFrameGenerator") << econd_header.size()
                                    << " word(s) of event packet header prepend. New size of ECON frame: "
                                    << econd_event.size();

    econd_event.push_back(computeCRC(econd_header));
    if (econd_params.add_idle_word) {
      const uint8_t buffer_status = 0, error_status = 0, reset_request = 0;
      econd_event.push_back(
          econd::buildIdleWord(buffer_status, error_status, reset_request, econd_params.programmable_pattern));
    }

    return econd_event;
  }

  std::vector<uint64_t> HGCalFrameGenerator::produceSlinkEvent(uint32_t fed_id) const {
    last_slink_emul_info_.clear();

    std::vector<uint64_t> slink_event;  // prepare the output S-link payload (will be "converted" to 128-bit)
    for (size_t i = 0; i < slink_.num_capture_blocks; ++i) {
      std::vector<uint64_t> cb_payload;
      std::vector<backend::ECONDPacketStatus> econd_statuses(max_num_econds_,
                                                             backend::ECONDPacketStatus::InactiveECOND);
      for (const auto& econd : econd_params_) {  // add all ECON-D payloads
        const auto& econd_id = econd.first;
        if (!econd.second.active) {
          econd_statuses[econd_id] = backend::ECONDPacketStatus::InactiveECOND;
          continue;
        }
        econd_statuses[econd_id] =
            backend::ECONDPacketStatus::Normal;  //TODO: also implement other ECON-D packet issues
        const auto econd_evt = to64bit(produceECONEvent(econd_id));
        cb_payload.insert(cb_payload.end(), econd_evt.begin(), econd_evt.end());
        last_slink_emul_info_.addECONDEmulatedInfo(econd_id, lastECONDEmulatedInfo());
      }
      const auto& eid = last_emul_event_.first;
      const uint64_t event_id = std::get<0>(eid), bx_id = std::get<1>(eid), orbit_id = std::get<2>(eid);
      const auto l1a_header = to64bit(backend::buildCaptureBlockHeader(bx_id, event_id, orbit_id, econd_statuses));
      LogDebug("HGCalFrameGenerator").log([&l1a_header](auto& log) { printWords(log, "l1a", l1a_header); });
      cb_payload.insert(cb_payload.begin(), l1a_header.begin(), l1a_header.end());  // prepend capture block header
      slink_event.insert(slink_event.end(),
                         cb_payload.begin(),
                         cb_payload.end());  // concatenate the capture block to the full S-link payload
    }

    if (slink_.store_header_trailer) {
      // build the S-link header words
      const uint32_t content_id = backend::buildSlinkContentId(backend::SlinkEmulationFlag::Subsystem, 0, 0);
      const auto& eid = last_emul_event_.first;
      const uint64_t event_id = std::get<0>(eid), bx_id = std::get<1>(eid), orbit_id = std::get<2>(eid);
      const auto slink_header =
          to128bit(backend::buildSlinkHeader(slink_.boe_marker, slink_.format_version, event_id, content_id, fed_id));
      slink_event.insert(slink_event.begin(), slink_header.begin(), slink_header.end());  // prepend S-link header

      // build the S-link trailer words
      const uint16_t daq_crc = 0, crc = 0;
      const uint32_t event_length = slink_event.size() - slink_header.size() - 1;
      const uint16_t status = backend::buildSlinkRocketStatus(false, false, false, false, false);
      const auto slink_trailer =
          to128bit(backend::buildSlinkTrailer(slink_.eoe_marker, daq_crc, event_length, bx_id, orbit_id, crc, status));
      slink_event.insert(slink_event.end(), slink_trailer.begin(), slink_trailer.end());  // append S-link trailer

      LogDebug("HGCalFrameGenerator").log([&slink_header, &slink_trailer](auto& log) {
        printWords(log, "slink header", slink_header);
        log << "\n";
        printWords(log, "slink trailer", slink_trailer);
      });
    }

    return slink_event;
  }

  uint32_t HGCalFrameGenerator::computeCRC(const std::vector<uint32_t>& event_header) const {
    uint32_t crc = 0x12345678;  //TODO: implement 32-bit Bluetooth CRC
    return crc;
  }
}  // namespace hgcal
