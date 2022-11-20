#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void export_pyddx_data(py::module& m) {
  py::module data =
        m.def_submodule("data", "Useful data (Unit conversion, atomic radii, ...)");

  // Units
  const double angstrom = 1 / 0.52917721092;  // conversion factor Angström to Bohr
  data.attr("angstrom") = angstrom;

  // UFF atomic radii
  // Atomic radii derived from the universal force field
  // A. K. Rappe et. al. J. Am. Chem. Soc., 1992, 114 (25), pp 10024-10035
  // https://doi.org/10.1021/ja00051a040, data given in Angström,
  // will be converted to Bohr. Note that keys are normalised to lower case.
  std::map<std::string, double> radius_uff{
        {"h", 1.4430},   //  1
        {"he", 1.8100},  //  2
        {"li", 1.2255},  //  3
        {"be", 1.3725},  //  4
        {"b", 2.0415},   //  5
        {"c", 1.9255},   //  6
        {"n", 1.8300},   //  7
        {"o", 1.7500},   //  8
        {"f", 1.6820},   //  9
        {"ne", 1.6215},  // 10
        {"na", 1.4915},  // 11
        {"mg", 1.5105},  // 12
        {"al", 2.2495},  // 13
        {"si", 2.1475},  // 14
        {"p", 2.0735},   // 15
        {"s", 2.0175},   // 16
        {"cl", 1.9735},  // 17
        {"ar", 1.9340},  // 18
        {"k", 1.9060},   // 19
        {"ca", 1.6995},  // 20
        {"sc", 1.6475},  // 21
        {"ti", 1.5875},  // 22
        {"v", 1.5720},   // 23
        {"cr", 1.5115},  // 24
        {"mn", 1.4805},  // 25
        {"fe", 1.4560},  // 26
        {"co", 1.4360},  // 27
        {"ni", 1.4170},  // 28
        {"cu", 1.7475},  // 29
        {"zn", 1.3815},  // 30
        {"ga", 2.1915},  // 31
        {"ge", 2.1400},  // 32
        {"as", 2.1150},  // 33
        {"se", 2.1025},  // 34
        {"br", 2.0945},  // 35
        {"kr", 2.0705},  // 36
        {"rb", 2.0570},  // 37
        {"sr", 1.8205},  // 38
        {"y", 1.6725},   // 39
        {"zr", 1.5620},  // 40
        {"nb", 1.5825},  // 41
        {"mo", 1.5260},  // 42
        {"tc", 1.4990},  // 43
        {"ru", 1.4815},  // 44
        {"rh", 1.4645},  // 45
        {"pd", 1.4495},  // 46
        {"ag", 1.5740},  // 47
        {"cd", 1.4240},  // 48
        {"in", 2.2315},  // 49
        {"sn", 2.1960},  // 50
        {"sb", 2.2100},  // 51
        {"te", 2.2350},  // 52
        {"i", 2.2500},   // 53
        {"xe", 2.2020},  // 54
        {"cs", 2.2585},  // 55
        {"ba", 1.8515},  // 56
        {"la", 1.7610},  // 57
        {"ce", 1.7780},  // 58
        {"pr", 1.8030},  // 59
        {"nd", 1.7875},  // 60
        {"pm", 1.7735},  // 61
        {"sm", 1.7600},  // 62
        {"eu", 1.7465},  // 63
        {"gd", 1.6840},  // 64
        {"tb", 1.7255},  // 65
        {"dy", 1.7140},  // 66
        {"ho", 1.7045},  // 67
        {"er", 1.6955},  // 68
        {"tm", 1.6870},  // 69
        {"yb", 1.6775},  // 70
        {"lu", 1.8200},  // 71
        {"hf", 1.5705},  // 72
        {"ta", 1.5850},  // 73
        {"w", 1.5345},   // 74
        {"re", 1.4770},  // 75
        {"os", 1.5600},  // 76
        {"ir", 1.4200},  // 77
        {"pt", 1.3770},  // 78
        {"au", 1.6465},  // 79
        {"hg", 1.3525},  // 80
        {"tl", 2.1735},  // 81
        {"pb", 2.1485},  // 82
        {"bi", 2.1850},  // 83
        {"po", 2.3545},  // 84
        {"at", 2.3750},  // 85
        {"rn", 2.3825},  // 86
        {"fr", 2.4500},  // 87
        {"ra", 1.8385},  // 88
        {"ac", 1.7390},  // 89
        {"th", 1.6980},  // 90
        {"pa", 1.7120},  // 91
        {"u", 1.6975},   // 92
        {"np", 1.7120},  // 93
        {"pu", 1.7120},  // 94
        {"am", 1.6905},  // 95
        {"cm", 1.6630},  // 96
        {"bk", 1.6695},  // 97
        {"cf", 1.6565},  // 98
        {"es", 1.6495},  // 99
        {"fm", 1.6430},  // 100
        {"md", 1.6370},  // 101
        {"no", 1.6240},  // 102
        {"lr", 1.6180},  // 103
  };
  for (auto& kv : radius_uff) kv.second *= angstrom;
  data.attr("radius_uff")         = radius_uff;  // Recommended radii are scale * uff
  data.attr("radius_uff_scaling") = 1.1;

  // Van der Waals radii taken from A. Bondi, J. Phys. Chem., 1964, 68, 441 and given
  // in Angström (converted to Bohr later). Keys are normalised to lower case.
  std::map<std::string, double> radius_bondi{
        {"h", 1.20},  {"he", 1.40}, {"li", 1.82}, {"c", 1.70},  {"n", 1.55},
        {"o", 1.52},  {"f", 1.47},  {"ne", 1.54}, {"na", 2.27}, {"mg", 1.73},
        {"si", 2.10}, {"p", 1.80},  {"s", 1.80},  {"cl", 1.75}, {"ar", 1.88},
        {"k", 2.75},  {"ni", 1.63}, {"cu", 1.40}, {"zn", 1.39}, {"ga", 1.87},
        {"as", 1.85}, {"se", 1.90}, {"br", 1.85}, {"kr", 2.02}, {"pd", 1.63},
        {"ag", 1.72}, {"cd", 1.58}, {"in", 1.93}, {"sn", 2.17}, {"te", 2.06},
        {"i", 1.98},  {"xe", 2.16}, {"pt", 1.75}, {"au", 1.66}, {"hg", 1.55},
        {"tl", 1.96}, {"pb", 2.02}, {"u", 1.86},
  };
  for (auto& kv : radius_bondi) kv.second *= angstrom;
  data.attr("radius_bondi") = radius_bondi;  // Recommended radii are scale * bondi
  data.attr("radius_bondi_scaling") = 1.2;

  // The full list of supported radii sets
  data.attr("radii_sets") = std::vector<std::string>{"uff", "bondi"};

  // Tabulated dielectric constants of solvents.
  // Taken from https://gaussian.com/scrf/. More can also be found at
  // https://comp.chem.umn.edu/solvation. Keys are normalised to lower case.
  data.attr("solvent_epsilon") = std::map<std::string, double>{
        {"water", 78.3553},
        {"acetonitrile", 35.688},
        {"methanol", 32.613},
        {"ethanol", 24.852},
        {"isoquinoline", 11.00},
        {"quinoline", 9.16},
        {"chloroform", 4.7113},
        {"diethylether", 4.2400},
        {"dichloromethane", 8.93},
        {"dichloroethane", 10.125},
        {"carbontetrachloride", 2.2280},
        {"benzene", 2.2706},
        {"toluene", 2.3741},
        {"chlorobenzene", 5.6968},
        {"nitromethane", 36.562},
        {"heptane", 1.9113},
        {"cyclohexane", 2.0165},
        {"aniline", 6.8882},
        {"acetone", 20.493},
        {"tetrahydrofuran", 7.4257},
        {"dimethylsulfoxide", 46.826},
        {"argon", 1.430},
        {"krypton", 1.519},
        {"xenon", 1.706},
        {"n-octanol", 9.8629},
        {"1,1,1-trichloroethane", 7.0826},
        {"1,1,2-trichloroethane", 7.1937},
        {"1,2,4-trimethylbenzene", 2.3653},
        {"1,2-dibromoethane", 4.9313},
        {"1,2-ethanediol", 40.245},
        {"1,4-dioxane", 2.2099},
        {"1-bromo-2-methylpropane", 7.7792},
        {"1-bromooctane", 5.0244},
        {"1-bromopentane", 6.269},
        {"1-bromopropane", 8.0496},
        {"1-butanol", 17.332},
        {"1-chlorohexane", 5.9491},
        {"1-chloropentane", 6.5022},
        {"1-chloropropane", 8.3548},
        {"1-decanol", 7.5305},
        {"1-fluorooctane", 3.89},
        {"1-heptanol", 11.321},
        {"1-hexanol", 12.51},
        {"1-hexene", 2.0717},
        {"1-hexyne", 2.615},
        {"1-iodobutane", 6.173},
        {"1-iodohexadecane", 3.5338},
        {"1-iodopentane", 5.6973},
        {"1-iodopropane", 6.9626},
        {"1-nitropropane", 23.73},
        {"1-nonanol", 8.5991},
        {"1-pentanol", 15.13},
        {"1-pentene", 1.9905},
        {"1-propanol", 20.524},
        {"2,2,2-trifluoroethanol", 26.726},
        {"2,2,4-trimethylpentane", 1.9358},
        {"2,4-dimethylpentane", 1.8939},
        {"2,4-dimethylpyridine", 9.4176},
        {"2,6-dimethylpyridine", 7.1735},
        {"2-bromopropane", 9.3610},
        {"2-butanol", 15.944},
        {"2-chlorobutane", 8.3930},
        {"2-heptanone", 11.658},
        {"2-hexanone", 14.136},
        {"2-methoxyethanol", 17.2},
        {"2-methyl-1-propanol", 16.777},
        {"2-methyl-2-propanol", 12.47},
        {"2-methylpentane", 1.89},
        {"2-methylpyridine", 9.9533},
        {"2-nitropropane", 25.654},
        {"2-octanone", 9.4678},
        {"2-pentanone", 15.200},
        {"2-propanol", 19.264},
        {"2-propen-1-ol", 19.011},
        {"3-methylpyridine", 11.645},
        {"3-pentanone", 16.78},
        {"4-heptanone", 12.257},
        {"4-methyl-2-pentanone", 12.887},
        {"4-methylpyridine", 11.957},
        {"5-nonanone", 10.6},
        {"aceticacid", 6.2528},
        {"acetophenone", 17.44},
        {"a-chlorotoluene", 6.7175},
        {"anisole", 4.2247},
        {"benzaldehyde", 18.220},
        {"benzonitrile", 25.592},
        {"benzylalcohol", 12.457},
        {"bromobenzene", 5.3954},
        {"bromoethane", 9.01},
        {"bromoform", 4.2488},
        {"butanal", 13.45},
        {"butanoicacid", 2.9931},
        {"butanone", 18.246},
        {"butanonitrile", 24.291},
        {"butylamine", 4.6178},
        {"butylethanoate", 4.9941},
        {"carbondisulfide", 2.6105},
        {"cis-1,2-dimethylcyclohexane", 2.06},
        {"cis-decalin", 2.2139},
        {"cyclohexanone", 15.619},
        {"cyclopentane", 1.9608},
        {"cyclopentanol", 16.989},
        {"cyclopentanone", 13.58},
        {"decalin-mixture", 2.196},
        {"dibromomethane", 7.2273},
        {"dibutylether", 3.0473},
        {"diethylamine", 3.5766},
        {"diethylsulfide", 5.723},
        {"diiodomethane", 5.32},
        {"diisopropylether", 3.38},
        {"dimethyldisulfide", 9.6},
        {"diphenylether", 3.73},
        {"dipropylamine", 2.9112},
        {"e-1,2-dichloroethene", 2.14},
        {"e-2-pentene", 2.051},
        {"ethanethiol", 6.667},
        {"ethylbenzene", 2.4339},
        {"ethylethanoate", 5.9867},
        {"ethylmethanoate", 8.3310},
        {"ethylphenylether", 4.1797},
        {"fluorobenzene", 5.42},
        {"formamide", 108.94},
        {"formicacid", 51.1},
        {"hexanoicacid", 2.6},
        {"iodobenzene", 4.5470},
        {"iodoethane", 7.6177},
        {"iodomethane", 6.8650},
        {"isopropylbenzene", 2.3712},
        {"m-cresol", 12.44},
        {"mesitylene", 2.2650},
        {"methylbenzoate", 6.7367},
        {"methylbutanoate", 5.5607},
        {"methylcyclohexane", 2.024},
        {"methylethanoate", 6.8615},
        {"methylmethanoate", 8.8377},
        {"methylpropanoate", 6.0777},
        {"m-xylene", 2.3478},
        {"n-butylbenzene", 2.36},
        {"n-decane", 1.9846},
        {"n-dodecane", 2.0060},
        {"n-hexadecane", 2.0402},
        {"n-hexane", 1.8819},
        {"nitrobenzene", 34.809},
        {"nitroethane", 28.29},
        {"n-methylaniline", 5.9600},
        {"n-methylformamide-mixture", 181.56},
        {"n,n-dimethylacetamide", 37.781},
        {"n,n-dimethylformamide", 37.219},
        {"n-nonane", 1.9605},
        {"n-octane", 1.9406},
        {"n-pentadecane", 2.0333},
        {"n-pentane", 1.8371},
        {"n-undecane", 1.9910},
        {"o-chlorotoluene", 4.6331},
        {"o-cresol", 6.76},
        {"o-dichlorobenzene", 9.9949},
        {"o-nitrotoluene", 25.669},
        {"o-xylene", 2.5454},
        {"pentanal", 10.0},
        {"pentanoicacid", 2.6924},
        {"pentylamine", 4.2010},
        {"pentylethanoate", 4.7297},
        {"perfluorobenzene", 2.029},
        {"p-isopropyltoluene", 2.2322},
        {"propanal", 18.5},
        {"propanoicacid", 3.44},
        {"propanonitrile", 29.324},
        {"propylamine", 4.9912},
        {"propylethanoate", 5.5205},
        {"p-xylene", 2.2705},
        {"pyridine", 12.978},
        {"sec-butylbenzene", 2.3446},
        {"tert-butylbenzene", 2.3447},
        {"tetrachloroethene", 2.268},
        {"tetrahydrothiophene-s,s-dioxide", 43.962},
        {"tetralin", 2.771},
        {"thiophene", 2.7270},
        {"thiophenol", 4.2728},
        {"trans-decalin", 2.1781},
        {"tributylphosphate", 8.1781},
        {"trichloroethene", 3.422},
        {"triethylamine", 2.3832},
        {"xylene-mixture", 2.3879},
        {"z-1,2-dichloroethene", 9.2},
  };
}
