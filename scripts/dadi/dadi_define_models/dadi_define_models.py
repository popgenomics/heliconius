#! ~/py33/bin
# -*- coding: utf-8 -*-

import numpy
import site
site.addsitedir('/data3/users/christec/dadi/dadi-1.6.3_modif')
import dadi



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"Defines the different models of divergence used in our inferences with DADI.

"There are six categories of models: "SI": Strict Isolation, "IM": Isolation with Migration, "AM": Ancient Migration, "PAM": Periodic Ancient Migration, "SC": Secondary Contact and "PSC": Periodic Secondary Contact. 

"Each model comes in subcategories: "no suffix": Reference model; "ex": Exponential growth; "2N": Background selection with 2 categories of population size in the genome (we assume population size is reduced by a given factor in non-recombining regions); "2M2P": Selection against migrants with 2 categories of migration rate in the genome (we assume migration is null in barrier regions); and the combinations between subcategories: "2N2M2P"; "2Nex"; "2M2Pex"; "2N2M2Pex".

"All models assume an initial ancestral population which is split in two daughter populations.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



""""""""""""""""""""""""""""""
"Strict Isolation models (SI)"
""""""""""""""""""""""""""""""

def SI(params, (n1,n2), pts):
    nu1, nu2, Ts, O = params
    """
    Model with split and strict isolation.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: Time of divergence in strict isolation.
    O: Proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # Ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Split event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Correctly oriented spectrum
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # Misoriented spectrum
    fsM = dadi.Numerics.reverse_array(fsO)

    # Final spectrum
    fs = O*fsO+(1-O)*fsM
    return fs


def SIex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, Ts, Te, O = params
    """
    Model with split and strict isolation; exponential growth.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after exponential growth.
    nu2: Size of population 2 after exponential growth.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in strict isolation.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    # Exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phi = dadi.Integration.two_pops(phi, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    # Correctly oriented spectrum
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    # Misoriented spectrum
    fsM = dadi.Numerics.reverse_array(fsO)

    # Final spectrum
    fs = O*fsO+(1-O)*fsM
    return fs


def SI2N(params, (n1,n2), pts):
    nu1, nu2, Ts, nr, bf, O = params
    """
    Model with split and strict isolation; two categories of population size in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ts: Time of divergence in strict isolation.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions. 
    O: Proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # Divergence in strict isolation
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
	# Correctly oriented spectrum
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
	# Misoriented spectrum
    fsnrM = dadi.Numerics.reverse_array(fsnrO) 
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    # Divergence in strict isolation
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
	# Correctly oriented spectrum
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
	# Misoriented spectrum
    fsrM = dadi.Numerics.reverse_array(fsrO)

    # Final spectrum
    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def SI2Nex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, Ts, Te, nr, bf, O = params
    """
    Model with split and strict isolation; exponential growth; two categories of population size in the genome. 

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after exponential growth.
    nu2: Size of population 2 after exponential growth.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in strict isolation.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    # Divergence in strict isolation
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    	# Exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinr = dadi.Integration.two_pops(phinr, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    	# Correctly oriented spectrum
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
	# Misoriented spectrum
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    # Divergence in strict isolation
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    	# Exponential growth
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phir = dadi.Integration.two_pops(phir, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    	# Correctly oriented spectrum
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    	# Misoriented spectrum
    fsrM = dadi.Numerics.reverse_array(fsrO)

    # Final spectrum
    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs



""""""""""""""""""""""""""""""""""""""
"Isolation with Migration models (IM)"
""""""""""""""""""""""""""""""""""""""

def IM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, O = params
    """
    Model with continuous migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in continuous migration.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Divergence in continuous migration
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def IMex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Te, O = params
    """ 
    Model with continuous migration; exponential growth. 
	
    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after exponential growth.
    nu2: Size of population 2 after exponential growth.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in strict isolation.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """ 
    xx = dadi.Numerics.default_grid(pts) 	
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=m12, m21=m21) 
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phi = dadi.Integration.two_pops(phi, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM 
    return fs


def IM2N(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, nr, bf, O = params
    """
    Model with continuous migration; two categories of population size in the genome. 

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.

    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=m12, m21=m21)
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def IM2Nex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Te, nr, bf, O = params
    """
    Model with continuous migration; exponential growth; two categories of population size in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after exponential growth.
    nu2: Size of population 2 after exponential growth.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in continuous migration.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinr = dadi.Integration.two_pops(phinr, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phir = dadi.Integration.two_pops(phir, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def IM2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, P1, P2, O = params
    """
    Model with continuous migration; two categories of migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in continuous migration.
    Te: Time of the exponential growth in continuous migration.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=m12, m21=0)
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=m21)
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs

    
def IM2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Te, P1, P2, O = params
    """
    Model with continuous migration; exponential growth; two categories of migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after exponential growth.
    nu2: Size of population 2 after exponential growth.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in continuous migration.
    Te: Time of the exponential growth in continuous migration.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def IM2N2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, nr, bf, P1, P2, O = params
    """
    Model with continuous migration; two categories of population size and migration rate in the genome. 

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=m12, m21=m21)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=m12, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=m21)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=m12, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=m21)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs


def IM2N2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Te, nr, bf, P1, P2, O = params
    """
    Model with continuous migration; exponential growth; two categories of population size and migration rate in the genome. 

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in continuous migration.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs



""""""""""""""""""""""""
"Ancient Migration (AM)"
""""""""""""""""""""""""

def AM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, O = params
    """
    Model with split, ancient migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Ancient migration event
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def AMex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, O = params
    """
    Model with split, ancient migration; exponential growth.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in strict isolation.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phi = dadi.Integration.two_pops(phi, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def AM2N(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, nr, bf, O = params
    """
    Model with split, ancient migration; two categories of population size in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def AM2Nex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, nr, bf, O = params
    """
    Model with split, ancient migration; exponential growth; two categories of population size in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinr = dadi.Integration.two_pops(phinr, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phir = dadi.Integration.two_pops(phir, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def AM2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, P1, P2, O = params
    """
    Model with split, ancient migration; two categories of migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def AM2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, P1, P2, O = params
    """
    Model with split, ancient migration; exponential growth; two categories of migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in continuous migration.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1a, nu2a, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)  
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1a, nu2a, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)  
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def AM2N2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, nr, bf, P1, P2, O = params
    """
    Model with split, ancient migration; two categories of population size and migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs


def AM2N2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, nr, bf, P1, P2, O = params
    """
    Model with split, ancient migration; exponential growth; two categories of population size and migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs



""""""""""""""""""""""""""""""""""
"Periodic Ancient Migration (PAM)"
""""""""""""""""""""""""""""""""""

def PAM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, O = params
    """
    Model with split, two periods of ancient migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    n1,n2: Size of fs to generate.
    O: The proportion of accurate SNP orientation.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Ancient migration event-1
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    # Ancient migration event-2
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    # Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def PAMex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, O = params
    """
    Model with split, two periods of ancient migration; exponential growth.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in strict isolation.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phi = dadi.Integration.two_pops(phi, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def PAM2N(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, nr, bf, O = params
    """
    Model with split, two periods of ancient migration; two categories of population size in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def PAM2Nex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, nr, bf, O = params
    """
    Model with split, two periods of ancient migration; exponential growth; two categories of population size in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tam, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinr = dadi.Integration.two_pops(phinr, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phir = dadi.Integration.two_pops(phir, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def PAM2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, P1, P2, O = params
    """
    Model with split, two periods of ancient migration; two categories of migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def PAM2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, P1, P2, O = params
    """
    Model with split, two periods of ancient migration; exponential growth; two categories of migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in continuous migration.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tam, nu1a, nu2a, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tam, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1a, nu2a, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tam, nu1a, nu2a, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)  
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1a, nu2a, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tam, nu1a, nu2a, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)  
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def PAM2N2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Tam, Ts, nr, bf, P1, P2, O = params
    """
    Model with split, two periods of ancient migration; two categories of population size and migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs


def PAM2N2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Tam, Ts, Te, nr, bf, P1, P2, O = params
    """
    Model with split, two periods of ancient migration; exponential growth; two categories of population size and migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Tam: Time of ancient migration.
    Ts: Time of divergence in strict isolation.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tam, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tam, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tam, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tam, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tam, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tam, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs



""""""""""""""""""""""""
"Secondary Contact (SC)"
""""""""""""""""""""""""

def SC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, O = params
    """
    Model with split, strict isolation, and secondary contact.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    #Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    #Secondary contact event
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def SCex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, O = params
    """
    Model with split, strict isolation, and secondary contact; exponential growth.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phi = dadi.Integration.two_pops(phi, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def SC2N(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, nr, bf, O = params
    """
    Model with split, strict isolation, and secondary contact; two categories of population size in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def SC2Nex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, nr, bf, O = params
    """
    Model with split, strict isolation, and secondary contact; exponential growth; two categories of population size in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinr = dadi.Integration.two_pops(phinr, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phir = dadi.Integration.two_pops(phir, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs= O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def SC2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, P1, P2, O = params
    """
    Model with split, strict isolation, and secondary contact; two categories of migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def SC2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, P1, P2, O = params
    """
    Model with split, strict isolation, and secondary contact; exponential growth; two categories of migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1a, nu2a, m12=0, m21=0)    
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1a, nu2a, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0) 
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1a, nu2a, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21) 
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def SC2N2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, nr, bf, P1, P2, O = params
    """
    Model with split, strict isolation, and secondary contact; two categories of population size and migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=m21)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs


def SC2N2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, nr, bf, P1, P2, O = params
    """
    Model with split, strict isolation, and secondary contact; exponential growth; two categories of population size and migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs



""""""""""""""""""""""""""""""""""
"Periodic Secondary Contact (PSC)"
""""""""""""""""""""""""""""""""""

def PSC(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    #Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    #Secondary contact event-1
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    #Divergence in strict isolation
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=0, m21=0)
    #Secondary contact event-2
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def PSCex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; exponential growth.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phi = dadi.Integration.two_pops(phi, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phi = dadi.Integration.two_pops(phi, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsO = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    fsM = dadi.Numerics.reverse_array(fsO)

    fs = O*fsO+(1-O)*fsM
    return fs


def PSC2N(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, nr, bf, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; two categories of population size in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1, nu2, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs = O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def PSC2Nex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, nr, bf, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; exponential growth; two categories of population size in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1.
    m21: Migration from population 1 to population 2.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
    phinr = dadi.PhiManip.phi_1D(xx)
    phinr = dadi.PhiManip.phi_1D_to_2D(xx, phinr)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    phinr = dadi.Integration.two_pops(phinr, xx, Ts, nu1a*bf, nu2a*bf, m12=0, m21=0)
    phinr = dadi.Integration.two_pops(phinr, xx, Tsc, nu1a*bf, nu2a*bf, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinr = dadi.Integration.two_pops(phinr, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsnrO = dadi.Spectrum.from_phi(phinr, (n1,n2), (xx,xx))
    fsnrM = dadi.Numerics.reverse_array(fsnrO)
    # Spectrum for recombining regions
    phir = dadi.PhiManip.phi_1D(xx)
    phir = dadi.PhiManip.phi_1D_to_2D(xx, phir)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    phir = dadi.Integration.two_pops(phir, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phir = dadi.Integration.two_pops(phir, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phir = dadi.Integration.two_pops(phir, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsrO = dadi.Spectrum.from_phi(phir, (n1,n2), (xx,xx))
    fsrM = dadi.Numerics.reverse_array(fsrO)

    fs= O*(nr*fsnrO + (1-nr)*fsrO) + (1-O) *(nr*fsnrM + (1-nr)*fsrM)
    return fs


def PSC2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, P1, P2, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; two categories of migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def PSC2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, P1, P2, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; exponential growth; two categories of migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-barrier regions in population 1 and 2
    phiN1N2 = dadi.PhiManip.phi_1D(xx)
    phiN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1N2)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Tsc, nu1a, nu2a, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1N2 = dadi.Integration.two_pops(phiN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21) 
    fsN1N2O = dadi.Spectrum.from_phi(phiN1N2, (n1,n2), (xx,xx))
    fsN1N2M = dadi.Numerics.reverse_array(fsN1N2O)
    # Spectrum for barrier regions in population 1 and 2
    phiI1I2 = dadi.PhiManip.phi_1D(xx)
    phiI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1I2)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1a, nu2a, m12=0, m21=0)    
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Tsc, nu1a, nu2a, m12=0, m21=0)    
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1I2 = dadi.Integration.two_pops(phiI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0) 
    fsI1I2O = dadi.Spectrum.from_phi(phiI1I2, (n1,n2), (xx,xx))
    fsI1I2M = dadi.Numerics.reverse_array(fsI1I2O)
    # Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phiN1I2 = dadi.PhiManip.phi_1D(xx)
    phiN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phiN1I2)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1a, nu2a, m12=m12, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Tsc, nu1a, nu2a, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiN1I2 = dadi.Integration.two_pops(phiN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0) 
    fsN1I2O = dadi.Spectrum.from_phi(phiN1I2, (n1,n2), (xx,xx))
    fsN1I2M = dadi.Numerics.reverse_array(fsN1I2O)
    # Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phiI1N2 = dadi.PhiManip.phi_1D(xx)
    phiI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phiI1N2)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1a, nu2a, m12=0, m21=m21)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Ts, nu1a, nu2a, m12=0, m21=0)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Tsc, nu1a, nu2a, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phiI1N2 = dadi.Integration.two_pops(phiI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21) 
    fsI1N2O = dadi.Spectrum.from_phi(phiI1N2, (n1,n2), (xx,xx))
    fsI1N2M = dadi.Numerics.reverse_array(fsI1N2O)

    fs = O*(P1*P2*fsN1N2O + (1-P1)*(1-P2)*fsI1I2O + P1*(1-P2)*fsN1I2O + (1-P1)*P2*fsI1N2O) + (1-O)*(P1*P2*fsN1N2M + (1-P1)*(1-P2)*fsI1I2M + P1*(1-P2)*fsN1I2M + (1-P1)*P2*fsI1N2M)
    return fs


def PSC2N2M2P(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts, Tsc, nr, bf, P1, P2, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; two categories of population size and migration rate in the genome.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=m21)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs


def PSC2N2M2Pex(params, (n1,n2), pts):
    nu1a, nu2a, nu1, nu2, m12, m21, Ts, Tsc, Te, nr, bf, P1, P2, O = params
    """
    Model with split, strict isolation, and two periods of secondary contact; exponential growth; two categories of population size and migration rate in the genome.

    nu1a: Size of population 1 after split.
    nu2a: Size of population 2 after split.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from population 2 to population 1 in non-barrier regions.
    m21: Migration from population 1 to population 2 in non-barrier regions.
    Ts: Time of divergence in strict isolation.
    Tsc: Time of secondary contact.
    Te: Time of the exponential growth in continuous migration.
    nr: Proportion of "non-recombining" regions affected by background selection.
    bf : Background factor, which defines the extent of population size reduction in "nr" regions.
    P1: Proportion of "non-barrier" regions in population 1.
    P2: Proportion of "non-barrier" regions in population 2.
    O: The proportion of accurate SNP orientation.
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    xx = dadi.Numerics.default_grid(pts)
    # Spectrum for non-recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phinrN1N2 = dadi.PhiManip.phi_1D(xx)
    phinrN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1N2)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1N2 = dadi.Integration.two_pops(phinrN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21)
    fsnrN1N2O = dadi.Spectrum.from_phi(phinrN1N2, (n1,n2), (xx,xx))
    fsnrN1N2M = dadi.Numerics.reverse_array(fsnrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phinrI1I2 = dadi.PhiManip.phi_1D(xx)
    phinrI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1I2)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1I2 = dadi.Integration.two_pops(phinrI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsnrI1I2O = dadi.Spectrum.from_phi(phinrI1I2, (n1,n2), (xx,xx))
    fsnrI1I2M = dadi.Numerics.reverse_array(fsnrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phinrN1I2 = dadi.PhiManip.phi_1D(xx)
    phinrN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrN1I2)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Tsc, nu1*bf, nu2*bf, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrN1I2 = dadi.Integration.two_pops(phinrN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsnrN1I2O = dadi.Spectrum.from_phi(phinrN1I2, (n1,n2), (xx,xx))
    fsnrN1I2M = dadi.Numerics.reverse_array(fsnrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phinrI1N2 = dadi.PhiManip.phi_1D(xx)
    phinrI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phinrI1N2)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=m21)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Ts, nu1*bf, nu2*bf, m12=0, m21=0)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Tsc, nu1*bf, nu2*bf, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phinrI1N2 = dadi.Integration.two_pops(phinrI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsnrI1N2O = dadi.Spectrum.from_phi(phinrI1N2, (n1,n2), (xx,xx))
    fsnrI1N2M = dadi.Numerics.reverse_array(fsnrI1N2O)
    # Spectrum for recombining regions
	# Spectrum for non-barrier regions in population 1 and 2
    phirN1N2 = dadi.PhiManip.phi_1D(xx)
    phirN1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1N2)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Tsc, nu1, nu2, m12=m12, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1N2 = dadi.Integration.two_pops(phirN1N2, xx, Te, nu1_func, nu2_func, m12=m12, m21=m21)
    fsrN1N2O = dadi.Spectrum.from_phi(phirN1N2, (n1,n2), (xx,xx))
    fsrN1N2M = dadi.Numerics.reverse_array(fsrN1N2O)
	# Spectrum for barrier regions in population 1 and 2
    phirI1I2 = dadi.PhiManip.phi_1D(xx)
    phirI1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1I2)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Tsc, nu1, nu2, m12=0, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1I2 = dadi.Integration.two_pops(phirI1I2, xx, Te, nu1_func, nu2_func, m12=0, m21=0)
    fsrI1I2O = dadi.Spectrum.from_phi(phirI1I2, (n1,n2), (xx,xx))
    fsrI1I2M = dadi.Numerics.reverse_array(fsrI1I2O)
	# Spectrum for non-barrier regions in population 1 and barrier regions in population 2
    phirN1I2 = dadi.PhiManip.phi_1D(xx)
    phirN1I2 = dadi.PhiManip.phi_1D_to_2D(xx, phirN1I2)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Tsc, nu1, nu2, m12=m12, m21=0)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirN1I2 = dadi.Integration.two_pops(phirN1I2, xx, Te, nu1_func, nu2_func, m12=m12, m21=0)
    fsrN1I2O = dadi.Spectrum.from_phi(phirN1I2, (n1,n2), (xx,xx))
    fsrN1I2M = dadi.Numerics.reverse_array(fsrN1I2O)
	# Spectrum for barrier regions in population 1 and non-barrier regions in population 2
    phirI1N2 = dadi.PhiManip.phi_1D(xx)
    phirI1N2 = dadi.PhiManip.phi_1D_to_2D(xx, phirI1N2)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Ts, nu1, nu2, m12=0, m21=0)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Tsc, nu1, nu2, m12=0, m21=m21)
    nu1_func = lambda t: numpy.exp(numpy.log(nu1) * t/Te)
    nu2_func = lambda t: numpy.exp(numpy.log(nu2) * t/Te)
    phirI1N2 = dadi.Integration.two_pops(phirI1N2, xx, Te, nu1_func, nu2_func, m12=0, m21=m21)
    fsrI1N2O = dadi.Spectrum.from_phi(phirI1N2, (n1,n2), (xx,xx))
    fsrI1N2M = dadi.Numerics.reverse_array(fsrI1N2O)

    fs = O*(nr*(P1*P2*fsnrN1N2O + (1-P1)*(1-P2)*fsnrI1I2O + P1*(1-P2)*fsnrN1I2O + (1-P1)*P2*fsnrI1N2O) + (1-nr)*(P1*P2*fsrN1N2O + (1-P1)*(1-P2)*fsrI1I2O + P1*(1-P2)*fsrN1I2O + (1-P1)*P2*fsrI1N2O)) + (1-O)*(nr*(P1*P2*fsnrN1N2M + (1-P1)*(1-P2)*fsnrI1I2M + P1*(1-P2)*fsnrN1I2M + (1-P1)*P2*fsnrI1N2M) + (1-nr)*(P1*P2*fsrN1N2M + (1-P1)*(1-P2)*fsrI1I2M + P1*(1-P2)*fsrN1I2M + (1-P1)*P2*fsrI1N2M))
    return fs

