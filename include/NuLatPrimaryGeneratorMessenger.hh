//	NuLatPrimaryGeneratorMessenger.hh 

#ifndef NuLatPrimaryGeneratorMessenger_h
#define NuLatPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class NuLatPrimaryParticleGenerator;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;

class NuLatPrimaryGeneratorMessenger: public G4UImessenger{
	public:

		NuLatPrimaryGeneratorMessenger(NuLatPrimaryParticleGenerator*);
		~NuLatPrimaryGeneratorMessenger();
    
	void SetNewValue(G4UIcommand*, G4String);
    
	private:

		NuLatPrimaryParticleGenerator*	NuLatGenerator;
		G4UIdirectory*					gunDir; 
		G4UIcmdWithAString* 			EventTypeCmd;
		G4UIcmdWithAString*				FileNameCmd;
		G4UIcmdWithAString*				InputCmd;
		G4UIcmdWithoutParameter*		UpdateCmd;
		std::string*					MessInput;
	
};


#endif

