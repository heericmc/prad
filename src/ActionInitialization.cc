#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new EventAction());
    SetUserAction(new StackingAction());   // kill secondary e-/e+/gamma for speed
}
