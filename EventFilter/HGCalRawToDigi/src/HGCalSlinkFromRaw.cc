#include "EventFilter/HGCalRawToDigi/interface/HGCalSlinkFromRaw.h"

// example reader by P.Dauncey, using https://gitlab.cern.ch/pdauncey/hgcal10glinkreceiver

using namespace hgcal;

SlinkFromRaw::SlinkFromRaw(const edm::ParameterSet &iConfig) : SlinkEmulatorBase(iConfig) {
 
  inputfiles_=iConfig.getUntrackedParameter<std::vector<std::string>>("inputs");
  ifile_=0;
    
  edm::LogInfo("SlinkFromRaw") << "files: \n";
  copy(begin(inputfiles_),end(inputfiles_), std::ostream_iterator<std::string>{std::cout,"\n"});
    
  // Make the buffer space for the records
  record_ = new hgcal_slinkfromraw::RecordT<4095>;
  nEvents_=0;
}

//
FEDRawDataCollection SlinkFromRaw::next() {

  FEDRawDataCollection raw_data;
  
  //open for the first time
  if( fileReader_.closed() ) {
    auto inputfile = inputfiles_[ifile_];
    fileReader_.open(inputfile);
    
  }

  //no more records in the file
  if(!fileReader_.read(record_)) {
    ifile_++;

    if(ifile_>=inputfiles_.size())
      throw cms::Exception("SlinkFromRaw::next") << "No more files";
    fileReader_.close();
    auto inputfile=inputfiles_[ifile_];
    fileReader_.open(inputfile);
  }

  edm::LogInfo("SlinkFromRaw: Reading record from file #") << ifile_ << "nevents=" << nEvents_ << "\n";
  
  // Set up specific records to interpet the formats
  const hgcal_slinkfromraw::RecordStarting *rStart((hgcal_slinkfromraw::RecordStarting*)record_);
  const hgcal_slinkfromraw::RecordStopping *rStop((hgcal_slinkfromraw::RecordStopping*)record_);
  const hgcal_slinkfromraw::RecordRunning  *rEvent((hgcal_slinkfromraw::RecordRunning*)record_);
  if(record_->state()==hgcal_slinkfromraw::FsmState::Starting) {
    rStart->print();
    std::cout << std::endl;
  } else if(record_->state()==hgcal_slinkfromraw::FsmState::Stopping){
    rStop->print();
    std::cout << std::endl;
  } else {                
    // We have a new event
    nEvents_++;
    bool print(nEvents_<=1);
    if(print) {
      rEvent->print();
      std::cout << std::endl;
    }

    // Check id is correct
    if(!rEvent->valid()) rEvent->print();
                
    // Access the Slink header ("begin-of-event")
    // This should always be present; check pattern is correct
    const hgcal_slinkfromraw::SlinkBoe *b(rEvent->slinkBoe());
    assert(b!=nullptr);
    if(!b->validPattern()) b->print();
    
    // Access the Slink trailer ("end-of-event")
    // This should always be present; check pattern is correct
    const hgcal_slinkfromraw::SlinkEoe *e(rEvent->slinkEoe());
    assert(e!=nullptr);
    if(!e->validPattern()) e->print();
                
    // Access the BE packet header
    const hgcal_slinkfromraw::BePacketHeader *bph(rEvent->bePacketHeader());
    if(bph!=nullptr && print) bph->print();
    
    // Access ECON-D packet as an array of 32-bit words
    const uint32_t *pEcond(rEvent->econdPayload());
    
    // Check this is not an empty event
    if(pEcond!=nullptr) {
                    
      if(print) {
        std::cout << "First 10 words of ECON-D packet" << std::endl;
        std::cout << std::hex << std::setfill('0');
        for(unsigned i(0);i<10;i++) {
          std::cout << "0x" << std::setw(8) << pEcond[i] << std::endl;
        }
        std::cout << std::dec << std::setfill(' ');
        std::cout << std::endl;
      }

      //record_->print();
      
      //copy to the event      
      auto *payload=record_->getPayload();
      auto payloadLength=record_->payloadLength()-2;
      std::cout << payload << " " << payloadLength/sizeof(uint32_t) << std::endl;
      size_t total_event_size = payloadLength/sizeof(char);
      std::cout << "\t -> " << total_event_size << std::endl;
      auto& fed_data = raw_data.FEDData(1);
      fed_data.resize(total_event_size);
      auto* ptr = fed_data.data();
      memcpy(ptr, (char*)payload, total_event_size);
    }
  }

  return raw_data;
}
