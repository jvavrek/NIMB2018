#ifndef eventActionMessenger_hh
#define eventActionMessenger_hh 1

#include "G4UImessenger.hh"

class eventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class eventActionMessenger : public G4UImessenger {
 public:
  eventActionMessenger(eventAction *EA);
  ~eventActionMessenger();

  void SetNewValue(G4UIcommand *, G4String);

 private:
  eventAction *EventAction;
  G4UIdirectory *detDataDirectory;
  G4UIcmdWithAnInteger *eventFreqCmd;
};

#endif
