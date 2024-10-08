import os
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")

from httcp.production.PolarimetricA1 import PolarimetricA1
#from IPython import embed


def GetPhiCP(
        p4hcandinfodict: dict, 
        method_leg1: str, 
        method_leg2: str,
        mode_leg1: str, 
        mode_leg2: str,
        **kwars
) -> ak.Array:
    #print(f"{method_leg1}-{method_leg2}")
    #print(f"{mode_leg1}-{mode_leg2}")
    """
    Inputs:
      p4hcandinfodict : full hcand dict
        e.g. 
          hcand : [h1 (i.e. e/mu/tau), h2 (i.e. tau)]
          p4hcandinfodict = {
            "p4h1"    : p4 of e/mu/tau i.e. the 1st component of hcand,
            "p4h1pi"  : p4 of charged pion if tau else empty,
            "p4h1pi0" : p4 of reconstructed pi0 if tau else empty,
            "p4h2"    : p4 of tau (always) i.e. the 2nd component of hcand,
            "p4h2pi"  : p4 of charged pion from tau,
            "p4h2pi0" : p4 of reconstructed pi0 from tau,
          }
      method_leg1     : "DP/PV/IP"
      method_leg2     : "DP/PV/IP"
      mode_leg1       : "e/mu/pi/rho/a1"
      mode_leg2       : "pi/rho/a1" because the 2nd component of hcand is always tauh
    Output:
      Final PhiCP array for any methods / configurations
    Steps:
      Prepares the input vectors for PhiCP
      Compute the PhiCP
    """
    #from IPython import embed; embed()
    phicp_input_vec_dict = PrepareVecsForPhiCP(p4hcandinfodict,
                                               method_leg1,
                                               method_leg2, 
                                               mode_leg1,
                                               mode_leg2)
    phicp = ak.values_astype(ComputeAcopAngle(phicp_input_vec_dict),
                             np.float32)
    phicp = ak.enforce_type(phicp, "var * float32")
    
    return phicp


def PrepareVecsForPhiCP(
        p4hcandinfodict: dict, 
        method_leg1: str, 
        method_leg2: str, 
        mode_leg1: str, 
        mode_leg2: str
) -> dict:
    """
    Inputs:
      p4hcandinfodict : full hcand dict
      method_leg1     : "DP/PV/IP"
      method_leg2     : "DP/PV/IP"
      mode_leg1       : "e/mu/pi/rho/a1"
      mode_leg2       : "pi/rho/a1"
    Output:
      DP method:
          P : Pion at zero momentum frame
          R : Pi0 (Transverse) at zero momentum frame
          H : R [same, it will be used to compute the sign]
          Y : Phase shift
          C : Charge of Tau or hcand1/2
      IP method:
          P : e/mu(for etau/mutau channel) or pion(for tautau channel) at zero momentum frame
          R : IP vector at zero momentum frame
          H : R
    Steps:
      Prepare input vectors for PhiCP calculation
    """
    p4h1    = p4hcandinfodict["p4h1"]
    p4h1pi  = p4hcandinfodict["p4h1pi"]
    p4h1pi0 = p4hcandinfodict["p4h1pi0"]
    p4h2    = p4hcandinfodict["p4h2"]
    p4h2pi  = p4hcandinfodict["p4h2pi"]
    p4h2pi0 = p4hcandinfodict["p4h2pi0"]
    
    _P1, _R1, _C1 = _prepareVecs(p4h1, p4h1pi, p4h1pi0, method_leg1, mode_leg1)
    _P2, _R2, _C2 = _prepareVecs(p4h2, p4h2pi, p4h2pi0, method_leg2, mode_leg2)

    # Get Boost
    boostv = _getBoost(_P1, _R1,
                       _P2, _R2,
                       method_leg1, 
                       method_leg2,
                       mode_leg1, 
                       mode_leg2)

    P1, R1, H1, Y1 = _reStructureVecs(boostv, _P1, _R1, p4h1pi0, method_leg1, mode_leg1)
    P2, R2, H2, Y2 = _reStructureVecs(boostv, _P2, _R2, p4h2pi0, method_leg2, mode_leg2)

    return {"P1" : P1, "R1" : R1, "H1": H1, "C1": _C1,
            "P2" : P2, "R2" : R2, "H2": H2, "C2": _C2,
            "Y1" : Y1, "Y2" : Y2}



def ComputeAcopAngle(vecsdict):
    """
    Geometrical estimation of PhiCP
    For PVPV:
       P1/P2  --> full tau p4
       R1/R2  --> k i.e. full tau x polarimetric vector
       H1/H2  --> h i.e. polarimetric vector
    For DPDP: 
       P1/P2  --> vecPiPlus
       R1/R2  --> vecPiZeroPlustransv
       H1/H2  --> R1/R2
       Y1     --> phase-shift for leg 1
       Y2     --> phase-shift for leg 2
    """
    #embed()
    assert len(vecsdict) == 10, "input dict to ComputeAcopAngle does not have proper structure"
    P1 = vecsdict["P1"]
    R1 = vecsdict["R1"]
    H1 = vecsdict["H1"]
    C1 = vecsdict["C1"]
    P2 = vecsdict["P2"]
    R2 = vecsdict["R2"]
    H2 = vecsdict["H2"]
    C2 = vecsdict["C2"]
    Y1 = vecsdict["Y1"]
    Y2 = vecsdict["Y2"]

    Y  = None
    if Y1 is not None:
        if Y2 is not None:
            Y = Y1*Y2
        else:
            Y = Y1
    else:
        if Y2 is not None:
            Y = Y2
            
    acop = np.arccos(R1.dot(R2))
    sign = _getSign(P1,H1,P2,H2,C1,C2) #P1.dot(H2.cross(H1))
    acop = ak.where(sign < 0.0, 2.*np.pi - acop, acop)
    if Y is not None:
        acop  = ak.where(Y < 0.0, acop + np.pi, acop)
        #acop  = ak.where((Y < 0.0) & (acop > 2.*np.pi), 
        #                 acop - 2.*np.pi, 
        #                 acop)

    #Map  angles into [0,2pi] interval
    acop = ak.where(acop > 2.*np.pi, acop - 2.* np.pi, acop) 
    acop = ak.where(acop < 0,        acop + 2.* np.pi, acop)
    return acop


# PRIVATE
def _getSign(P1,H1,P2,H2,C1,C2):
    #embed()
    Pm = ak.where(C1 < 0, P1, P2) # check tau-
    Hm = ak.where(C1 < 0, H1, H2) # Sort according to tau charge
    Hp = ak.where(C1 < 0, H2, H1) # Same
    sign = Pm.dot(Hp.cross(Hm))
    #sign = P2.dot(H1.cross(H2))
    return sign
    

def _getBoost(
        P1   : ak.Array, 
        R1   : ak.Array, 
        P2   : ak.Array, 
        R2   : ak.Array, 
        leg1_method: str, 
        leg2_method: str, 
        leg1_mode: str, 
        leg2_mode: str,
        **kwargs
) -> ak.Array:
    """
    A function to get the boostvec of the corresponding frame
    To get the frame:
      for leg1_method == "DP" and leg2_method == "DP":
        frame = P1 + P2 i.e. Pion from tau1 and Pion from tau2
      for leg1_method == "IP" and leg2_method == "IP":  
        frame = P1 + P2 i.e. e/mu/pion from tau1 and pion from tau2

    returns:
      boostvec of the frame
    """
    frame = None
    if leg1_method == "PV" and leg2_method == "PV":
        # P1/P2 -- tau1/2, R1/R2 -- pion_tau1/2
        if leg1_mode == "pi" and leg2_mode == "pi":
            frame = R1 + R2
        else:
            frame = P1 + P2
    else:
        frame = P1 + P2

    return frame.boostvec


def _prepareVecs(
        hcand: ak.Array, 
        hcand_pi: ak.Array, 
        hcand_pi0: ak.Array, 
        leg_method: str, 
        leg_mode: str
) -> tuple[ak.Array, ak.Array, ak.Array]:
    """
    A private function to be used in PrepareVecsForPhiCP
    This mainly returns the required vectors for different modes and methods
    Returned vectors will further be modified in the reStructure step
    -- Output:
    Method : PV
      Mode : Pi, Rho, a1 [only possible modes]
        WARNING: Frame will be different for Pi
                 and it will use R instead of P
        _P --> tau (hcand)
        _R --> pi  (pions)
        It will be converted to polarimertric vector in next function
    Method : DP
      Mode : Rho, a1 [only possible modes]
        _P --> pi  (diff for a1)
        _R --> pi0 (diff for a1)
    Method : IP
      Mode : e, mu, pi [only needed]
        _P --> pion if mode == pi else e/mu
        _R --> e/mu/tau
      
    """
    _P = None
    _R = None

    if leg_method == "DP":
        # https://github.com/alebihan/IPHCProductionTools/blob/master/HiggsCPinTauDecays/TauDecaysInterface/src/SCalculator.cc#L363-L366
        if leg_mode == "rho":
            _P = hcand_pi      # charged pion
            _R = hcand_pi0     # pi0
        elif leg_mode == "a1":
            _P = _get_pi_a1_DP(hcand_pi) # check invMass of (osPi + ss1Pi) and (osPi + ss2Pi), then select the one closest to rho mass
            _R = hcand_pi[:,:1]          # osPion
        else:
            raise RuntimeError(f"Wrong mode : {leg_mode}")


    elif leg_method == "PV":
        _P = hcand
        _R = hcand_pi

    elif leg_method == "IP":
        if leg_mode == "e" or leg_mode == "mu":
            _P = hcand
            _R = hcand
        elif leg_mode == "pi":
            _P = hcand_pi
            _R = hcand
        #pass

    else:
        raise RuntimeError(f"Wrong {leg_method}")

    #embed()
    hcand_charge = hcand.charge
    #print(f"leg_method: {leg_method}, leg_mode: {leg_mode}")
    #print(f"P: {P}")
    #print(f"R: {R}")
    return _P, _R, hcand_charge


def _reStructureVecs(
        boostv: ak.Array, 
        V1: ak.Array, 
        V2: ak.Array, 
        V3: ak.Array, 
        leg_method: str, 
        leg_mode: str
) -> tuple[ak.Array, ak.Array, ak.Array, ak.Array]:
    """
    PostProcess the inputs to phicp with boost
    Output : 
    Method : PV
      Mode : Pi, Rho, a1 [only possible modes]
        WARNING: Frame will be different for Pi
        P --> tau (hcand)
        R --> k   (tau x polarimetric vector)
        H --> h   (pv)
    Method : DP
      Mode : Rho, a1 [only possible modes]
        P --> charged pion at ZMF  (diff for a1: checking the closest one by requiring rho invm)
        R --> transverse component of pi0 at ZMF (diff for a1: get the os pi)
        H --> pi0 (same as R)
    Method : IP
      Mode : e, mu, pi [only needed]
        P --> e/mu/pion at ZMF
        R --> transverse component of IP at ZMF
        H --> R
    """
    P = None
    R = None
    H = None # to calculate angle
    Y = None
    
    
    if leg_method == "DP":
        """
        Input
          V1   : Pi   [_P from _prepareVecs]
          V2   : Pi0  [_R from _prepareVecs]
          V3   : hcand_pi0 [not required]
        Output
          P : Pion at zero momentum frame
          R : Pi0 (Transverse) at zero momentum frame
          H : R [same, as it will be used to compute the sign]
          Y : Phase shift
        """
        pi_ZMF  = V1.boost(boostv.negative())
        #embed()

        pi0_ZMF = V2.boost(boostv.negative())
        pi_ZMF_unit  = pi_ZMF.pvec.unit
        pi0_ZMF_unit = pi0_ZMF.pvec.unit
        pi0_ZMF_unit_T = (pi0_ZMF_unit - pi_ZMF_unit*(pi_ZMF_unit.dot(pi0_ZMF_unit))).unit
    
        P = pi_ZMF_unit
        R = pi0_ZMF_unit_T
        H = R
        Y = (pi_ZMF.energy - pi0_ZMF.energy)/(pi_ZMF.energy + pi0_ZMF.energy)
    
    elif leg_method == "PV":
        """
        Input
          V1   : tau [i.e. hcand]
          V2   : Pi  [hcand decay]
          V3   : Pi0 [hcand decay]
        """

        P, R, H = _pv(boostv, V1, V2, V3, leg_mode)

    elif leg_method == "IP":
        """
        V1   : p4 of e/mu/pi
        V2   : p4 of e/mu/tau
        V3   : p4 of Pi0 [not relevant]
        """
        P_ZMF = V1.boost(boostv.negative()).pvec
        # build the IP vector and then apply negative boost
        IP = ak.zip({"x":V2.IPx, "y":V2.IPy, "z":V2.IPz, "t":ak.zeros_like(V2.IPz)},
                    with_name="LorentzVector",
                    behavior=coffea.nanoevents.methods.vector.behavior)
        R_ZMF      = IP.boost(boostv.negative()).pvec

        P_ZMF_unit   = P_ZMF.unit
        R_ZMF_unit   = R_ZMF.unit
        #R_ZMF_T_unit = R_ZMF_unit.add((P_ZMF_unit.multiply(R_ZMF_unit.dot(P_ZMF_unit))).negative())
        R_ZMF_T_unit = (R_ZMF_unit - P_ZMF_unit*(R_ZMF_unit.dot(P_ZMF_unit))).unit
        
        #P_ZMF_unit = P_ZMF.unit
        #R_ZMF_T = (R_ZMF - (R_ZMF.dot(P_ZMF))*P_ZMF_unit)
        #R_ZMF_T_unit = (R_ZMF - (R_ZMF.dot(P_ZMF))*P_ZMF_unit).unit
        
        P = P_ZMF_unit
        R = R_ZMF_T_unit
        H = R
        #pass


    else:
        raise RuntimeError(f"Wrong method: {leg_method}")

    return P, R, H, Y


def _get_pi_a1_DP(p4_pi):
    Minv1 = (p4_pi[:,:1] + p4_pi[:,1:2]).mass
    Minv2 = (p4_pi[:,:1] + p4_pi[:,2:3]).mass
    Pi = ak.where(np.abs(0.77526-Minv1) < np.abs(0.77526-Minv2), p4_pi[:,1:2], p4_pi[:,2:3])
    return Pi


def _pv(boostv: ak.Array, 
        tau : ak.Array, 
        pi  : ak.Array, 
        pi0 : ak.Array, 
        leg_mode: str) -> ak.Array:
    P  = tau.boost(boostv.negative())
    if leg_mode == "pi":
        pv = pi.boost(boostv.negative()).pvec
    elif leg_mode == "rho":
        pi  = pi.boost(boostv.negative())
        pi0 = pi0.boost(boostv.negative())
        q   = pi.subtract(pi0)
        N   = P.subtract(pi.add(pi0))
        pv  = (((2*(q.dot(N))*q.pvec).subtract(q.mass2*N.pvec)))
    elif leg_mode == "a1":
        os_pi_HRF   = pi[:, 0:1].boost(boostv.negative())
        ss1_pi_HRF  = pi[:, 1:2].boost(boostv.negative()) 
        ss2_pi_HRF  = pi[:, 2:3].boost(boostv.negative()) 
        a1pol       = PolarimetricA1(P,
                                     os_pi_HRF,
                                     ss1_pi_HRF,
                                     ss2_pi_HRF,
                                     tau.charge)
        pv          = -a1pol.PVC().pvec
    else:
        raise RuntimeError(f"Wrong mode: {leg_mode}")
    
    P = P.pvec.unit
    H = pv.unit
    R = (H.cross(P)).unit

    return P, R, H
