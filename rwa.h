/******************************************************************************************/
/*  rwa.h v_1.0 (MIS_based rwa solution) - Seyed Jafar Mosavian, Ruijia 9/11/2017         */
/******************************************************************************************/
#include "head.h"
#include "rwa_struct.h"

#ifdef _MSC_VER
#define DLL_EXPORT __declspec( dllexport ) 
#else
#define DLL_EXPORT
#endif

//number is bigger when the print information is more important.
//#define PRINT_DETAIL_0
//#define PRINT_DETAIL_1
#define PRINT_DETAIL_2

//output path&wave assignment information stored in pathList.
DLL_EXPORT void rwa(RWA_DB *db);

//config[path_num]: output configuration;
//dual[1+reqGroupNum]: input dual values;
//link_list: array of pointers that point to a set of links;
//path_list: array of pointers that point to a set of paths;
//return: cost.
double pricing_problem_path(RWA_DB *db, double config[], const double dual[], const int link_num, const int path_num,
                              RWA_LINK *link_list[], RWA_PATH *path_list[]);

//config: point to config_list[num_configurations];
//return: 0-if z_c is binary and its wavelength is used by all links. 
int column_generation(RWA_DB *db, const double dual[], double *cost, WAVE_CONFIG *config);

//config_list[]: input configurations;
//num_configurations: input current number of existing configurations;
//dual[1+reqGroupNum]: output dualvalues;
//return: sum(y_sd).
int master_problem(RWA_DB *db, WAVE_CONFIG config_list[], const int num_configurations, const int W_used, double dual[]);

//transfer solution into path&lambda results.
void process_solutions(RWA_DB *db, const WAVE_CONFIG config_list[], const int num_configurations, const int used_wave_num);

//return: number of binary z_c.
int preset_config_list(RWA_DB *db, WAVE_CONFIG config_list[]);

void test_instance_1();
void test_instance_2();
void test_instance_3();