/******************************************************************************************/
/*  rwa_struct.h v_1.0 (MIS_based rwa solution) - Seyed Jafar Mosavian, Ruijia 9/11/2017  */
/******************************************************************************************/
typedef struct wave_config {
    double z_c;//chosen number of each configuration;
    int is_binary;//1- z_c could only be 1/0; 0- z_c could be integer.
    int lambda_index;//1<=index<=W, valid only if is_binary=1.
    double *config;//config array, len=pathNum.
    int *ac_sd;//ac_sd[reqGroupNum].
} WAVE_CONFIG;

typedef struct rwa_link{
    int src;
    int des;
    int index;
    int *wave;//1-used; 0-free.
} RWA_LINK;

void init_link(RWA_LINK *link, int index, int src, int des, int W);
void print_link_list(RWA_LINK *linkList, int linkNum);
void del_link_list(RWA_LINK *linkList, int linkNum);
void set_link_wave(RWA_LINK *link, int wave_index, int wave_status);
void print_link_wave(RWA_LINK *linkList, int linkNum, int W);

typedef struct rwa_path{
    int index;
    int src;
    int des;
    int link_num;
    int req_index;
    int *link_list;
    int *wave;//1-used; 0-free.

    //calc result.
    int assigned_wave_num;
    int *assigned_wave;//1 if assigned wave, len = W.
} RWA_PATH;

void init_path(RWA_PATH *path, int index, int src, int des, int link_num, int link_list[], int W);
void print_path_list(RWA_PATH *pathList, int pathNum);
void del_path_list(RWA_PATH *pathList, int pathNum);
void set_path_wave(RWA_PATH pathList[], int pathNum, RWA_LINK linkList[], int W);
void reset_path_result(RWA_PATH pathList[], int pathNum, int W);
void print_path_result(RWA_PATH pathList[], int pathNum, int W);
void print_path_wave(RWA_PATH *pathList, int pathNum, int W);

typedef struct rwa_req{
    int index;
    int src;
    int des;
    int req_num;
    int path_num;
    int *path_list;
} RWA_REQ;

void init_req(RWA_REQ *req, int index, int src, int des, int req_num, int path_num, int path_list[]);
void print_req_list(RWA_REQ *reqList, int reqGroupNum);
void del_req_list(RWA_REQ *reqList, int reqGroupNum);
void set_path_req_index(RWA_REQ *reqList, int reqGroupNum, RWA_PATH *pathList, int pathNum);

typedef struct rwa_db {
    int W;
    int linkNum;
    int pathNum;
    int reqGroupNum;
    RWA_LINK *linkList;//link db.
    RWA_PATH *pathList;//path db.
    RWA_REQ *reqList;//requests in sd pair db.
} RWA_DB;

void init_rwa_db(RWA_DB *db, int linkNum, int pathNum, int reqGroupNum, int W);
void del_db(RWA_DB *db);
