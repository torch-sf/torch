"""
Torch parameter file implementation

Joshua Wall, Drexel University
To be re-organized by Aaron Tran, aaron.tran@columbia.edu

Implementation TBD . . . do we want

(1) plaintext par file like FLASH,
(2) stand-alone python file w/ hooks to AMUSe, so user can do interesting ICs,
(3) something else?

In any case, we need a scheme to handle defaults cleanly.
Ex. if we do a python file, we should probably do a simple class inheritance...

possibly...

    class TorchPar(TorchDefaultPar):
        def __init__(self, blah):
            pass

or, we could do...

    from torch import main()
    from torch_par import TorchPar()

    par = TorchPar()  # inherits from abstract base calss

    par['min_sf_mass'] = 1.0 | units.Msun
    par['nproc'] = 42  # user can put their own submit script parser here
    par['tmax'] = 1e13

    run = TorchRun(par)  # user might interface with "run" as well.
    run.init()
    run.go()


    # user could swap out a different par file after a pre-determined time.
    par['tmax'] = 2e13
    par['use_radiation'] = True

    run.set_par(par)
    run.go()

    run.clean()

    # user can query state, print information, dump diagnostics

in principle state could be configured other ways... have to think about this.

beware making too complex.  this is not possible within a summer, needs at
least a month or two of dedicated work...
"""




if __name__ == '__main__':
    pass
