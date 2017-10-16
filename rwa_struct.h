/******************************************************************************************/
/*  rwa_struct.h v_3.0 (MIS_based rwa solution) - Seyed Jafar Mosavian, Ruijia 10/8/2017  */
/******************************************************************************************/
typedef struct wave_config{
    double z_c;//chosen number of each configuration;
    int is_binary;//1- z_c could only be 1/0; 0- z_c could be integer.
    int lambda_index;//1<=index<=W, valid only if is_binary=1.
    double *config;//config array, len=pathNum.
    int *ac_sd;//ac_sd[reqGroupNum].
} WAVE_CONFIG;

typedef struct rwa_link{
    int index;
	int src;
    int des;
    int *assigned_wavelength;// 0-free, 1-used.
	int *lightpath_list;//deployed lightpaths
} RWA_LINK;

void init_link(RWA_LINK *link, int index, int src, int des, int W);
void print_link_list(RWA_LINK *linkList, int linkNum);
void del_link_list(RWA_LINK *linkList, int linkNum);
void set_link_wavelength(RWA_LINK *link, int wave_index, int wave_status);
void print_link_wavelenth(RWA_LINK *linkList, int linkNum, int W);

typedef struct rwa_path{
    int index;
    int src;
    int des;
    int link_num;
	int *link_list;
    int sub_req_index;
	int *available_wavelength;// 0-free, 1-busy.
//    int usage_num;// 0-free, integer-used.
} RWA_PATH;

void init_path(RWA_PATH *path, int index, int src, int des, int link_num, int link_list[], int W);
void print_path_list(RWA_PATH *pathList, int pathNum);
void del_path_list(RWA_PATH *pathList, int pathNum);
void set_path_available_wavelength(RWA_PATH pathList[], int pathNum, RWA_LINK linkList[], int W);
void print_path_wavelength(RWA_PATH *pathList, int pathNum, int W);

typedef struct rwa_lightpath{
    int index;
    int sub_path_num;// predefined = 1, it will be incremented if there is any realy node(s)through the route.
    RWA_PATH *sub_path_list;
    int *assigned_wavelength_index;
	//int req_index;
} RWA_LIGHTPATH;

void init_lightpath(RWA_LIGHTPATH *lightpath, int index, int subPathNum, RWA_PATH subPathList[], int assignedWavelengthIndex[]);
void print_lightpath_list(RWA_LIGHTPATH *lightpathList, int lightpathNum);
void del_lightpath_list(RWA_LIGHTPATH *lightpathList, int lightpathNum);


typedef struct rwa_req{
    int index;
    int src;
    int des;
    int req_num;
    int path_num;
    int *path_list;// list of index of paths
} RWA_REQ;

void init_req(RWA_REQ *req, int index, int src, int des, int req_num, int path_num, int path_list[]);
void print_req_list(RWA_REQ *reqList, int reqGroupNum);
void del_req_list(RWA_REQ *reqList, int reqGroupNum);
void set_path_req_index(RWA_REQ *reqList, int reqGroupNum, RWA_PATH *pathList, int pathNum);

typedef struct rwa_db {
    int W;
    int link_num;
    int path_num;
	int lightpath_num;
    int req_qroup_num;
	int lightpath_index;// number of solved lightpaths
    RWA_LINK *link_list;//link db.
    RWA_PATH *path_list;//path db.
    RWA_REQ *req_list;//requests in sd pair db.
	RWA_LIGHTPATH *lightpath_list;//lightpath db.
} RWA_DB;

void init_rwa_db(RWA_DB *db, int linkNum, int pathNum, int lightpathNum, int reqGroupNum, int W);
void reset_rwa_db(RWA_DB *db, int pathNum, int reqGroupNum);
void del_db(RWA_DB *db);

void reset_lightpath_links(RWA_DB *db, RWA_LIGHTPATH *lightpath);
