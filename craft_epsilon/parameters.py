"""The 41 fitted parameters of CRAFT-epsilon (ALPB diethyl ether conformer ensembles)."""



# Fitted in three stages to 1393 liquids (ALPB diethyl ether ensembles): nonpolar (3), single- and
# multiple-polar (21), hydrogen-bonding (17). Atom groups without a fitted scale keep 1; H is fixed at 1.
PARAMETERS = {
    "nonpolar": {                                   # stage 1: 3 parameters
        "atom_group_scales": {"C_SP2": 1.136800937817024, "C_SP3": 1.1194088762590442},
        "dipole_scale": 0.5163655332828004,         # c_np: mu_eff = c_np * mu
    },
    "polar": {                                      # stage 2: 21 parameters (atom scales add to stage 1)
        "atom_group_scales": {
            "Br": 0.3739473897481668, "C_SP": 1.0817948653402307, "Cl": 2.2910990264826854,
            "F": 1.727368827193598, "N_SP": 5.0065700381810805, "N_SP2": 4.325755097559717,
            "O_SP2": 1.0314509378419374, "O_SP3": 1.1414530077413967, "P_SP3": 1.6828188833753703,
            "S_SP3": 0.8866465009918179, "Si_SP3": 1.0240760404384397},
        "family_scales": {
            "carbonyl": 0.7917879423730423, "halogen": 0.5377559693139289, "nitro": 0.5722955213185961,
            "nitrogen": 0.7317787451883259, "other": 0.9447354609727306, "oxygen": 0.9491047229900367,
            "phosphorus": 0.9778962583087751, "silicon": 0.23383563546773975, "sulfur": 0.9407888410039368},
        "rest_scale": 0.43706219970398746,          # b
    },
    "hbond": {                                      # stage 3: 17 parameters (atom scales of stage 2)
        "fragment_scales": {
            "acetamide": 1.5900916124735547, "alcohol": 1.2593857616558493, "amide": 1.7012738852254643,
            "amine": 0.748395620969103, "carbonate/carbamate": 0.7126603876495889,
            "carboxylic_acid": 0.3113913566182205, "chloro": 1.4534340499851775e-13,
            "ester": 0.7626391086369851, "ether": 4.9983118741458516e-18, "fluoro": 0.7755049380975545,
            "formamide": 1.2554184560017758, "ketone": 0.9349410810306614, "nitro": 0.7331499569668267,
            "phenol": 0.8694745285854373, "thiol": 0.6643926820428278},
        "rare_scale": 0.9430214787799336,           # fragment types in < 2 training molecules and unseen types
        "rest_scale": 2.0124979804581414e-14,       # b_HB
    },
}
