#ifndef rootStorageMessenger_hh
#define rootStorageMessenger_hh 1

#include "G4UImessenger.hh"

class rootStorageManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;


class rootStorageManagerMessenger : public G4UImessenger {
 public:
  rootStorageManagerMessenger(rootStorageManager *);
  ~rootStorageManagerMessenger();

  void SetNewValue(G4UIcommand *, G4String);

 private:
  rootStorageManager *theRSManager;

  G4UIdirectory *ZKDirectory, *rootDirectory;

  G4UIcmdWithAString *rootFileNameCmd;
  G4UIcmdWithoutParameter *rootInitCmd;
  G4UIcmdWithoutParameter *rootWriteCmd;
};

#endif
