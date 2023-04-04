#include <cstddef>

#include "EventFilter/HGCalRawToDigi/interface/HGCalRawDataPackingTools.h"
#include "EventFilter/HGCalRawToDigi/interface/HGCalRawDataDefinitions.h"
#include "EventFilter/HGCalRawToDigi/interface/SlinkTypes.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

std::vector<uint32_t> hgcal::econd::produceERxData(
    const std::vector<bool>& chmap, const ERxData& erx, bool passZS, bool passZSm1, bool hasToA, bool char_mode) {
  auto format_word = [&erx, &passZS, &passZSm1, &hasToA, &char_mode](size_t i) -> std::pair<uint32_t, uint8_t> {
    if (i >= erx.tctp.size())
      throw cms::Exception("HGCalEmulator")
          << "Not enough channels in eRx data payload: " << i << " >= " << erx.tctp.size() << ".";
    if (char_mode)  // characterisation mode is easy
      return std::make_pair(((uint8_t)erx.tctp.at(i) & 0x3) << 30 | (erx.adc.at(i) & 0x3ff) << 20 |
                                (erx.tot.at(i) & 0x3ff) << 10 | (erx.toa.at(i) & 0x3ff),
                            32);
    switch (erx.tctp.at(i)) {
      case ToTStatus::ZeroSuppressed: {
        if (!passZS)
          throw cms::Exception("HGCalEmulator") << "ToT status is ZeroSuppressed, but event frame does not pass ZS.";
        if (passZSm1) {
          if (hasToA)
            return std::make_pair((0x1 << 30) | ((erx.adcm.at(i) & 0x3ff) << 20) | ((erx.adc.at(i) & 0x3ff) << 10) |
                                      (erx.toa.at(i) & 0x3ff),
                                  32);
          else
            return std::make_pair(((erx.adcm.at(i) & 0x3ff) << 10) | (erx.adc.at(i) & 0x3ff), 24);
        }
        // at this point, does not have any BX-1 ZS info
        if (hasToA)
          return std::make_pair((0x3 << 20) | ((erx.adc.at(i) & 0x3ff) << 10) | (erx.toa.at(i) & 0x3ff), 24);
        return std::make_pair((0x1 << 10) | (erx.adc.at(i) & 0x3ff), 16);
      }
      case ToTStatus::noZeroSuppressed_TOASuppressed:
        return std::make_pair((0x2 << 20) | ((erx.adcm.at(i) & 0x3ff) << 10) | (erx.adc.at(i) & 0x3ff), 24);
      case ToTStatus::invalid:
        return std::make_pair(
            (0x2 << 30) | ((erx.adcm.at(i) & 0x3ff) << 20) | ((erx.adc.at(i) & 0x3ff) << 10) | (erx.toa.at(i) & 0x3ff),
            32);
      case ToTStatus::AutomaticFull:
        return std::make_pair(
            (0x3 << 30) | ((erx.adcm.at(i) & 0x3ff) << 20) | ((erx.tot.at(i) & 0x3ff) << 10) | (erx.toa.at(i) & 0x3ff),
            32);
      default:
        throw cms::Exception("HGCalEmulator")
            << "Invalid ToT status retrieved for channel " << i << ": " << (int)erx.tctp.at(i) << ".";
    }
  };

  std::vector<uint32_t> data{0};
  std::vector<uint32_t>::iterator it_data = data.begin();

  uint8_t msb = 0;
  for (size_t i = 0; i < chmap.size(); ++i) {
    if (!chmap.at(i))
      continue;
    auto [word, nbits] = format_word(i);  // retrieve the channel word in its proper formatting

    if (msb >= 32)  // if we are already at the next word
      it_data = data.insert(data.end(), 0);
    msb %= 32;  // 0 <= msb < 32

    if (msb > 0)  // do we have some room for additional information?
      *it_data &= ((1 << msb) - 1);
    if (msb + nbits > 32) {  // spilling onto the next word
      uint8_t nbits_word1 = 32 - msb - nbits;
      *it_data |= (word & ((1 << nbits_word1) - 1)) << msb;
      it_data = data.insert(data.end(), word >> nbits_word1);
    } else  // everything fits into one word
      *it_data |= word << msb;
    msb += nbits;
  }
  return data;
}

std::vector<uint32_t> hgcal::econd::addChannelData(uint8_t& msb,
                                                   uint8_t tctp,
                                                   uint16_t adc,
                                                   uint16_t tot,
                                                   uint16_t adcm,
                                                   uint16_t toa,
                                                   bool passZS,
                                                   bool passZSm1,
                                                   bool hasToA,
                                                   bool charmode) {
  msb %= 32;

  if (charmode)  // characterization mode is easy
    return {uint32_t((tctp & 0x3) << 30 | (adc & 0x3ff) << 20 | (tot & 0x3ff) << 10 | (toa & 0x3ff))};

  std::vector<uint32_t> chan_words;
  // normal mode will depend on the TcTp and ZS states
  uint32_t raw_word = 0;
  uint8_t nbits = 0;

  switch (tctp) {
    case 0x0:
      if (passZS && passZSm1 && !hasToA) {
        raw_word = ((adcm & 0x3ff) << 10) | (adc & 0x3ff);
        nbits = 24;
      } else if (passZS && !passZSm1 && !hasToA) {
        raw_word = (1 << 10) | (adc & 0x3ff);
        nbits = 16;
      } else if (passZS && !passZSm1 && hasToA) {
        raw_word = (3 << 20) | ((adc & 0x3ff) << 10) | (toa & 0x3ff);
        nbits = 24;
      } else if (passZS && passZSm1 && hasToA) {
        raw_word = (1 << 30) | ((adcm & 0x3ff) << 20) | ((adc & 0x3ff) << 10) | (toa & 0x3ff);
        nbits = 32;
      }
      break;
    case 0x1:
      raw_word = (2 << 20) | ((adcm & 0x3ff) << 10) | (adc & 0x3ff);
      nbits = 24;
      break;
    case 0x2:
      raw_word = (2 << 30) | ((adcm & 0x3ff) << 20) | ((adc & 0x3ff) << 10) | (toa & 0x3ff);
      nbits = 32;
      break;
    case 0x3:
      raw_word = (3 << 30) | ((adcm & 0x3ff) << 20) | ((tot & 0x3ff) << 10) | (toa & 0x3ff);
      nbits = 32;
      break;
    default:
      throw cms::Exception("HGCalEmulator") << "Invalid TcTp mode: " << tctp << ".";
  }

  // pack raw word to 32b starting at the required msb
  // if it still fits, it sits in the current word
  if (nbits + msb <= 32) {
    uint32_t mask = 0xffffffff >> (32 - nbits);
    if (nbits > 0) {
      chan_words.push_back((raw_word & mask) << (32 - msb - nbits));
      msb += nbits;
    }
    return chan_words;
  }

  //if it doesn't fit start a new word
  //    32                         0 32                            0 32
  //    |           .    msb        |  nbits-msb  . 32-(nbits-msb)  |
  uint8_t spill(nbits - msb);
  uint32_t rawWord1(raw_word >> spill);
  chan_words.push_back(rawWord1);

  uint32_t rawWord2(raw_word << (32 - spill));
  chan_words.push_back(rawWord2);

  msb = spill != 32 ? 32 - spill : 32;
  return chan_words;
}

//
std::vector<uint32_t> hgcal::econd::eRxSubPacketHeader(
    uint8_t stat, uint8_t ham, bool bitE, uint16_t cm0, uint16_t cm1, std::vector<bool> chmap) {
  uint64_t chmap64b(0);
  size_t i = 0;
  for (const auto& ch : chmap) {
    chmap64b |= (ch << i);
    ++i;
  }
  return hgcal::econd::eRxSubPacketHeader(stat, ham, bitE, cm0, cm1, chmap64b);
}

//
std::vector<uint32_t> hgcal::econd::eRxSubPacketHeader(
    uint8_t stat, uint8_t ham, bool bitE, uint16_t cm0, uint16_t cm1, uint64_t chmap) {
  std::vector<uint32_t> header = {(stat & hgcal::ECOND_FRAME::ERXSTAT_MASK) << hgcal::ECOND_FRAME::ERXSTAT_POS |
                                  (ham & hgcal::ECOND_FRAME::ERXHAM_MASK) << hgcal::ECOND_FRAME::ERXHAM_POS |
                                  (cm0 & hgcal::ECOND_FRAME::COMMONMODE0_MASK) << hgcal::ECOND_FRAME::COMMONMODE0_POS |
                                  (cm1 & hgcal::ECOND_FRAME::COMMONMODE1_MASK) << hgcal::ECOND_FRAME::COMMONMODE1_POS};

  //summarize the channel status map
  uint32_t chmapw0(chmap & hgcal::ECOND_FRAME::CHMAP0_MASK);
  uint32_t chmapw1((chmap >> 32) & hgcal::ECOND_FRAME::CHMAP32_MASK);

  //add the channel map
  if (chmapw0 == 0 && chmapw1 == 0) {  // empty channels map (empty eRx)
    header[0] |= (bitE << hgcal::ECOND_FRAME::ERX_E_POS);
    header[0] |= (1 << hgcal::ECOND_FRAME::ERXFORMAT_POS);  //raise the F bit (empty eRX)
  } else {
    header[0] |= (chmapw1 & hgcal::ECOND_FRAME::CHMAP32_MASK) << hgcal::ECOND_FRAME::CHMAP32_POS;
    header.push_back((chmapw0 & hgcal::ECOND_FRAME::CHMAP0_MASK) << hgcal::ECOND_FRAME::CHMAP0_POS);
  }

  return header;
}

static uint8_t computeCRC8bit(const std::vector<uint32_t>& event_header) {
  uint8_t crc = 0x42;  //TODO: implement 8-bit Bluetooth CRC
  return crc;
}

//
std::vector<uint32_t> hgcal::econd::eventPacketHeader(uint16_t header,
                                                      uint16_t payload,
                                                      bool bitP,
                                                      bool bitE,
                                                      uint8_t ht,
                                                      uint8_t ebo,
                                                      bool bitM,
                                                      bool bitT,
                                                      uint8_t ehHam,
                                                      uint16_t bx,
                                                      uint16_t l1a,
                                                      uint8_t orb,
                                                      bool bitS,
                                                      uint8_t rr) {
  std::vector<uint32_t> words(2, 0);

  words[0] = (header & hgcal::ECOND_FRAME::HEADER_MASK) << hgcal::ECOND_FRAME::HEADER_POS |
             (payload & hgcal::ECOND_FRAME::PAYLOAD_MASK) << hgcal::ECOND_FRAME::PAYLOAD_POS |
             bitP << hgcal::ECOND_FRAME::BITP_POS | bitE << hgcal::ECOND_FRAME::BITE_POS |
             (ht & hgcal::ECOND_FRAME::HT_MASK) << hgcal::ECOND_FRAME::HT_POS |
             (ebo & hgcal::ECOND_FRAME::EBO_MASK) << hgcal::ECOND_FRAME::EBO_POS |
             bitM << hgcal::ECOND_FRAME::BITM_POS | bitT << hgcal::ECOND_FRAME::BITT_POS |
             (ehHam & hgcal::ECOND_FRAME::EHHAM_MASK) << hgcal::ECOND_FRAME::EHHAM_POS;

  words[1] = (bx & hgcal::ECOND_FRAME::BX_MASK) << hgcal::ECOND_FRAME::BX_POS |
             (l1a & hgcal::ECOND_FRAME::L1A_MASK) << hgcal::ECOND_FRAME::L1A_POS |
             (orb & hgcal::ECOND_FRAME::ORBIT_MASK) << hgcal::ECOND_FRAME::ORBIT_POS |
             bitS << hgcal::ECOND_FRAME::BITS_POS | (rr & hgcal::ECOND_FRAME::RR_MASK) << hgcal::ECOND_FRAME::RR_POS;

  const auto crc = computeCRC8bit(words);
  words[1] |= (crc & hgcal::ECOND_FRAME::EHCRC_MASK) << hgcal::ECOND_FRAME::EHCRC_POS;

  return words;
}

uint32_t hgcal::econd::buildIdleWord(uint8_t bufStat, uint8_t err, uint8_t rr, uint32_t progPattern) {
  return (progPattern & 0xffffff) << 8 | (rr & 0x3) << 6 | (err & 0x7) << 3 | (bufStat & 0x7) << 0;
}

//
std::vector<uint32_t> hgcal::backend::buildCaptureBlockHeader(
    uint32_t bc, uint32_t ec, uint32_t oc, const std::vector<hgcal::backend::ECONDPacketStatus>& econd_statuses) {
  if (econd_statuses.size() > 12)
    throw cms::Exception("HGCalEmulator") << "Invalid size for ECON-D statuses: " << econd_statuses.size() << ".";
  std::vector<uint32_t> header(2, 0);
  header[0] = (bc & hgcal::BACKEND_FRAME::CAPTUREBLOCK_BC_MASK) << hgcal::BACKEND_FRAME::CAPTUREBLOCK_BC_POS |
              (ec & hgcal::BACKEND_FRAME::CAPTUREBLOCK_EC_MASK) << hgcal::BACKEND_FRAME::CAPTUREBLOCK_EC_POS |
              (oc & hgcal::BACKEND_FRAME::CAPTUREBLOCK_OC_MASK) << hgcal::BACKEND_FRAME::CAPTUREBLOCK_OC_POS |
              (econd_statuses[11] & 0x7) << 1 | ((econd_statuses[10] >> 2) & 0x1);
  for (size_t i = 0; i < 11; i++)
    header[1] |= (econd_statuses[i] & 0x7) << i * 3;
  return header;
}

//
std::vector<uint32_t> hgcal::backend::buildSlinkHeader(
    uint8_t boe, uint8_t v, uint64_t global_event_id, uint32_t content_id, uint32_t fed_id) {
  std::vector<uint32_t> header(4, 0);
  header[0] = (boe & hgcal::BACKEND_FRAME::SLINK_BOE_MASK) << hgcal::BACKEND_FRAME::SLINK_BOE_POS |
              (v & hgcal::BACKEND_FRAME::SLINK_V_MASK) << hgcal::BACKEND_FRAME::SLINK_V_POS |
              ((global_event_id >> 41) & SLINK_GLOBAL_EVENTID_MSB_MASK)
                  << hgcal::BACKEND_FRAME::SLINK_GLOBAL_EVENTID_MSB_POS;
  header[1] = (global_event_id & SLINK_GLOBAL_EVENTID_LSB_MASK);
  header[2] = (content_id & hgcal::BACKEND_FRAME::SLINK_CONTENTID_MASK) << hgcal::BACKEND_FRAME::SLINK_CONTENTID_POS;
  header[3] = (fed_id & hgcal::BACKEND_FRAME::SLINK_SOURCEID_MASK) << hgcal::BACKEND_FRAME::SLINK_SOURCEID_POS;

  return header;
}

//
std::vector<uint32_t> hgcal::backend::buildSlinkTrailer(uint8_t eoe,
                                                        uint16_t daqcrc,
                                                        uint32_t event_length,
                                                        uint16_t bxid,
                                                        uint32_t orbit_id,
                                                        uint16_t crc,
                                                        uint16_t status) {
  std::vector<uint32_t> trailer(4, 0);

  trailer[0] = (eoe & hgcal::BACKEND_FRAME::SLINK_EOE_MASK) << hgcal::BACKEND_FRAME::SLINK_EOE_POS |
               (daqcrc & hgcal::BACKEND_FRAME::SLINK_DAQCRC_MASK) << hgcal::BACKEND_FRAME::SLINK_DAQCRC_POS;
  trailer[1] = (event_length & hgcal::BACKEND_FRAME::SLINK_EVLENGTH_MASK) << hgcal::BACKEND_FRAME::SLINK_EVLENGTH_POS |
               (bxid & hgcal::BACKEND_FRAME::SLINK_BXID_MASK) << hgcal::BACKEND_FRAME::SLINK_BXID_POS;
  trailer[2] = (orbit_id & hgcal::BACKEND_FRAME::SLINK_ORBID_MASK) << hgcal::BACKEND_FRAME::SLINK_ORBID_POS;
  trailer[3] = (crc & hgcal::BACKEND_FRAME::SLINK_CRC_MASK) << hgcal::BACKEND_FRAME::SLINK_CRC_POS |
               (status & hgcal::BACKEND_FRAME::SLINK_STATUS_MASK) << hgcal::BACKEND_FRAME::SLINK_STATUS_POS;

  return trailer;
}

//
uint32_t hgcal::backend::buildSlinkContentId(SlinkEmulationFlag e, uint8_t l1a_subtype, uint16_t l1a_fragment_cnt) {
  return 0x0 | (l1a_fragment_cnt & 0xffff) | (l1a_subtype & 0xff) << 16 | (e & 0x3) << 24;
}

uint16_t hgcal::backend::buildSlinkRocketStatus(
    bool fed_crc_err, bool slinkrocket_crc_err, bool source_id_err, bool sync_lost, bool fragment_trunc) {
  return 0x0 | (fed_crc_err & 0x1) << 0 | (slinkrocket_crc_err & 0x1) << 1 | (source_id_err & 0x1) << 2 |
         (sync_lost & 0x1) << 3 | (fragment_trunc & 0x1) << 4;
}
