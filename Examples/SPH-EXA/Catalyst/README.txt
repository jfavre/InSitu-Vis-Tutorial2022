The source code file catalyst_adaptor.h defines a single channel called "grid". Thus, all Python code will get the registered object of the same name,via the call:

reader = TrivialProducer(registrationName="grid")

a SLURM script is provided as example to execute the simulation and the associated ParaView Catalyst Python script
