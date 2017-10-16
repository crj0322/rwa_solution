#include "dynamic_rwa.h"

void d_rwa(RWA_DB *db)
{
    //CAUTION!!! every array/matrix index starts from 0, except for the glpk arrays: ia[] ja[] ar[].
    int ii, jj;
    int num_configurations = 0;
    double cost;
    double *dual;
    WAVE_CONFIG *config_list;
    int used_wave_num;
    int sum_ysd = 0, sum_req = 0;
    int ret = 0;

    set_path_wave(db->pathList, db->pathNum, db->linkList, db->W);
#ifdef PRINT_DETAIL_2
    print_req_list(db->reqList, db->reqGroupNum);
    print_path_list(db->pathList, db->pathNum);
    print_link_list(db->linkList, db->linkNum);
    print_path_wave(db->pathList, db->pathNum, db->W);
    print_path_result(db->pathList, db->pathNum, db->W);
    //print_link_wave(db->linkList, db->linkNum, db->W);
#endif

    for (ii=0; ii<db->reqGroupNum; ii++)
    {
        sum_req += db->reqList[ii].req_num;
    }

    dual = (double *)malloc((db->reqGroupNum+1) * sizeof(double));
    config_list = (WAVE_CONFIG *)malloc(db->W * sizeof(WAVE_CONFIG));

    used_wave_num = d_preset_config_list(db, config_list);

    //init dual values.
    dual[0] = 0;
    for (ii=0; ii<db->reqGroupNum; ii++)
    {
        dual[ii+1] = 1;
    }

    //solve.
    do 
    {
        if (num_configurations == db->W)
        {
            break;
        }
        //generate wavelength configuration.
        config_list[num_configurations].config = (double *)malloc(db->pathNum * sizeof(double));
        for (ii=0; ii<db->pathNum; ii++)
        {
            config_list[num_configurations].config[ii] = 0;
        }
        ret = d_column_generation(db, dual, &cost, &config_list[num_configurations]);
        if (ret == 0)
        {
            num_configurations++;
            continue;
        }
        
        if (cost <= 0)
        {
            break;
        }

        //calc ac_sd[reqGroupNum].
        config_list[num_configurations].ac_sd = (int *)malloc(db->reqGroupNum * sizeof(int));
        for (ii=0; ii<db->reqGroupNum; ii++)
        {
            config_list[num_configurations].ac_sd[ii] = 0;
        }
        for (ii=0; ii<db->pathNum; ii++)
        {
            if (config_list[num_configurations].config[ii])
            {
                config_list[num_configurations].ac_sd[db->pathList[ii].req_index - 1]++;
            }
        }
#ifdef PRINT_DETAIL_0
        //print ac_sd.
        printf("ac_sd:\n");
        for (ii=0; ii<db->reqGroupNum; ii++)
        {
            for (jj=0; jj<=num_configurations; jj++)
            {
                printf("%d\t", config_list[jj].ac_sd[ii]);
            }
            printf("\n");
        }
#endif

        sum_ysd = d_master_problem(db, config_list, ++num_configurations, used_wave_num, dual);

    } while ((cost > 0) && (sum_ysd < sum_req));

    d_process_solutions(db, config_list, num_configurations, used_wave_num);

    //free configs.
    for (ii=0; ii<db->W;ii++)
    {
        if (config_list[ii].ac_sd != NULL)
        {
            free(config_list[ii].ac_sd);
            config_list[ii].ac_sd = NULL;
        }
        if (config_list[ii].config != NULL)
        {
            free(config_list[ii].config);
            config_list[ii].config = NULL;
        }
    }
    free(config_list);
    free(dual);

#ifdef PRINT_DETAIL_2
    print_path_result(db->pathList, db->pathNum, db->W);
    print_path_wave(db->pathList, db->pathNum, db->W);
    print_link_wave(db->linkList, db->linkNum, db->W);
#endif

    return;
}

int d_preset_config_list(RWA_DB *db, WAVE_CONFIG config_list[])
{
    int ii, jj, kk;

    for (ii=0; ii<db->W; ii++)
    {
        config_list[ii].is_binary = 0;
        config_list[ii].lambda_index = 0;
        config_list[ii].z_c = 0;
        config_list[ii].config = NULL;
        config_list[ii].ac_sd = NULL;
    }

    kk=0;//index of config_list.
    for (ii=0; ii<db->W; ii++)
    {
        //if any link had used this wavelength, add a binary z_c.
        for (jj=0; jj<db->linkNum; jj++)
        {
            if (db->linkList[jj].wave[ii])
            {
                config_list[kk].is_binary = 1;
                config_list[kk++].lambda_index = ii+1;
                break;
            }
        }
    }
    return kk;
}

void d_process_solutions(RWA_DB *db, const WAVE_CONFIG config_list[], const int num_configurations, const int used_wave_num)
{
    int ii, jj, kk;
    int start_wave_index;//start_wave_index for every configuration.
    RWA_LINK *link;
    int *wave_status;//1-used by any link, 0-free for all links;
    int *wave_map;//index=index in config_list, value=db wave_index-1.

#ifdef PRINT_DETAIL_2
    //print solution.
    printf("solution:\n");
    printf("z_c\t");
    for (jj=0; jj<num_configurations; jj++)
    {
        printf("%d\t", (int)config_list[jj].z_c);
    }
    printf("\n\n");
    for (ii=0; ii< db->pathNum; ii++)
    {
        printf("path[%d] ", ii+1);
        for (jj=0; jj<num_configurations; jj++)
        {
            printf("%d\t", (int)config_list[jj].config[ii]);
        }
        printf("\n");
    }
    printf("\n");
#endif

    wave_status = (int *)malloc(db->W * sizeof(int));
    wave_map = (int *)malloc(db->W * sizeof(int));

    //calc wave_map.
    for (ii=0; ii<db->W; ii++)
    {
        wave_status[ii] = 0;
    }
    for (ii=0; ii<used_wave_num; ii++)
    {
        wave_map[ii] = config_list[ii].lambda_index-1;
        wave_status[wave_map[ii]] = 1;
    }
    for (ii=0, jj=used_wave_num; ii<db->W; ii++)
    {
        if (wave_status[ii] == 0)
        {
            wave_map[jj++] = ii;
        }
    }

    printf("wave_map:");
    for (ii=0; ii<db->W; ii++)
    {
        printf(" %d", wave_map[ii]);
    }
    printf("\n");

    if (jj != ii)
    {
        printf("calc wave_map error!\n");
    }

    for (ii=0; ii<db->pathNum; ii++)
    {
        for (jj=0; jj<used_wave_num; jj++)
        {
            if ((config_list[jj].z_c) && (config_list[jj].config[ii]))
            {
                db->pathList[ii].assigned_wave_num++;
                db->pathList[ii].assigned_wave[config_list[jj].lambda_index-1] = 1;
            }
        }
        start_wave_index = used_wave_num;//index of wave_map.
        for (jj=used_wave_num; jj<num_configurations; jj++)
        {
            if ((config_list[jj].z_c) && (config_list[jj].config[ii]))
            {
                for (kk=0; kk<(int)config_list[jj].z_c; kk++)
                {
                    db->pathList[ii].assigned_wave_num++;
                    db->pathList[ii].assigned_wave[wave_map[start_wave_index+kk]] = 1;
                }
            }
            if (config_list[jj].z_c)
            {
                start_wave_index += (int)(config_list[jj].z_c);
            }
        }

        //set wave status to linkList.
        for (jj=0; jj<db->pathList[ii].link_num; jj++)
        {
            link = &db->linkList[db->pathList[ii].link_list[jj]-1];
            if (db->pathList[ii].assigned_wave_num > 0)
            {
                for (kk=0; kk<db->W; kk++)
                {
                    if (link->wave[kk] && db->pathList[ii].assigned_wave[kk])
                    {
                        printf("path%d link%d wave%d conflict!!!\n", ii+1, db->pathList[ii].link_list[jj], kk+1);
                        //process conflict.
                        //to do...
                    }
                    link->wave[kk] |= db->pathList[ii].assigned_wave[kk];
                }
            }
        }
    }

    free(wave_status);
    free(wave_map);

    set_path_wave(db->pathList, db->pathNum, db->linkList, db->W);

    return;
}

int d_master_problem(RWA_DB *db, WAVE_CONFIG config_list[], const int num_configurations, const int W_used, double dual[])
{
    int ii, jj;
    glp_prob *mp;
    const int rows = 1 + db->reqGroupNum, columns = num_configurations + db->reqGroupNum;
    double obj;
    int sum_ysd = 0;

    /* CAUTION!!! Every array data starts at element [1] instead of [0] - the element [0] belived to have its column or row's information. */
    int ia[1+LP_SIZE] = {0}, ja[1+LP_SIZE] = {0};
    double ar[1+LP_SIZE] = {0};

    mp = glp_create_prob();

    glp_set_prob_name(mp, "(Restricted) Master Problem");
    glp_set_obj_dir(mp, GLP_MAX);
    glp_add_rows(mp, rows);

    glp_set_row_name(mp, 1, "sum_c(z_c)<=W");
    glp_set_row_bnds(mp, 1, GLP_UP, 0.0, db->W - W_used);

    for (ii = 2; ii < 2 + db->reqGroupNum; ii++)
    {
        glp_set_row_name(mp, ii, "y[sd]-sum_c(ac[sd]*z_c)<=0");
        glp_set_row_bnds(mp, ii, GLP_UP, 0.0, 0.0);
    }

    glp_add_cols(mp, columns);
    for (ii = 1; ii < 1 + db->reqGroupNum; ii++)
    {
        glp_set_col_name(mp, ii, "y_sd");
        glp_set_col_bnds(mp, ii, GLP_DB, 0.0, db->reqList[ii - 1].req_num);
        glp_set_col_kind(mp, ii, GLP_IV);
        glp_set_obj_coef(mp, ii, 1.0);
    }

    for (ii = 1 + db->reqGroupNum; ii < 1 + columns; ii++)//num_configurations
    {
        glp_set_col_name(mp, ii, "z_c");
        if (config_list[ii-db->reqGroupNum-1].is_binary)
        {
            glp_set_col_kind(mp, ii, GLP_BV);
        }
        else
        {
            glp_set_col_bnds(mp, ii, GLP_LO, 0.0, 0.0);
            glp_set_col_kind(mp, ii, GLP_IV);
        }
        glp_set_obj_coef(mp, ii, 0.0);
    }

    ii = 1;
    for (jj = 1; jj < 1 + db->reqGroupNum; jj++)
    {
        ia[jj + columns * (ii - 1)] = ii;
        ja[jj + columns * (ii - 1)] = jj;
        ar[jj + columns * (ii - 1)] = 0; // a[i,j] = 0 for y_sd.
    }
    for (jj = 1 + db->reqGroupNum; jj < 1 + columns; jj++)//num_configurations
    {
        ia[jj + columns * (ii - 1)] = ii;
        ja[jj + columns * (ii - 1)] = jj;
        ar[jj + columns * (ii - 1)] = 1; // a[i,j] = 1 for z_c.
    }

    for (ii = 2; ii < 2 + db->reqGroupNum; ii++)
    {
        for (jj = 1; jj < 1 + db->reqGroupNum; jj++)
        {
            ia[jj + columns * (ii - 1)] = ii;
            ja[jj + columns * (ii - 1)] = jj;
            if (jj == ii - 1)
            {
                ar[jj + columns * (ii - 1)] = 1; // a[i,j] = 1 for y_sd[jj].
            }
        }
        for (jj = 1 + db->reqGroupNum; jj < 1 + columns; jj++)//num_configurations
        {
            ia[jj + columns * (ii - 1)] = ii;
            ja[jj + columns * (ii - 1)] = jj;
            ar[jj + columns * (ii - 1)] = -(double)config_list[jj - db->reqGroupNum - 1].ac_sd[ii - 2]; // a[i,j] = -ac_sd[i,j].
        }
    }

#ifdef PRINT_DETAIL_0
    printf("row bounds:\n");
    printf("sum_c(z_c)<=W\n");
    for (ii = 2; ii < 2 + db->reqGroupNum; ii++)
    {
        printf("y_sd[%d]-sum_c(ac[%d]*z_c)<=0\n", ii - 1, ii - 1);
    }
    printf("var bounds:\n");
    for (ii = 1; ii < 1 + db->reqGroupNum; ii++)
    {
        printf("0<=y_sd[%d]<=%d\n", ii, db->reqList[ii - 1].req_num);
    }
    for (ii = 1 + db->reqGroupNum; ii < 1 + columns; ii++)//num_configurations
    {
        printf("z_c[%d]>=0\n", ii - db->reqGroupNum);
    }
    printf("s.t. coef:\n");
    for (ii=1; ii<1+rows; ii++)
    {
        for (jj=1; jj<1+columns; jj++)
        {
            printf("%.1f\t", ar[jj + columns * (ii - 1)]);
        }
        printf("\n");
    }
#endif

    glp_load_matrix(mp, rows * columns, ia, ja, ar);
    glp_simplex(mp, NULL);
    glp_intopt(mp, NULL);
    obj = glp_mip_obj_val(mp);

    for (ii = 1; ii < 1 + db->reqGroupNum; ii++)
    {
        sum_ysd += (int)glp_mip_col_val(mp, ii);
    }
    for (ii = 1 + db->reqGroupNum; ii < 1 + columns; ii++)//num_configurations
    {
        config_list[ii - db->reqGroupNum - 1].z_c = glp_mip_col_val(mp, ii);
    }

    //get dual values.
    dual[0] = glp_get_row_dual(mp, 1);
    for (ii = 2; ii < 2 + db->reqGroupNum; ii++)
    {
        dual[ii-1] = glp_get_row_dual(mp, ii);
    }

#ifdef PRINT_DETAIL_1
    printf("\nMASTER obj = %g", obj);
    for (ii = 1; ii < 1 + db->reqGroupNum; ii++)
    {
        printf("\n y[%d] = %g", ii, glp_mip_col_val(mp, ii));
    }
    printf("\n");
    for (ii = 1 + db->reqGroupNum; ii < 1 + columns; ii++)//num_configurations
    {
        printf("\n z_c[%d] = %g", ii - db->reqGroupNum, config_list[ii - db->reqGroupNum - 1].z_c);
    }
    printf("\n");
#endif

    return sum_ysd;
}

int d_column_generation(RWA_DB *db, const double dual[], double *cost, WAVE_CONFIG *config)
{
    int link_num = 0, path_num = 0, ii, jj;
    RWA_LINK **tmp_link_list;
    RWA_PATH **tmp_path_list;
    double *tmp_config;
    int lambda = config->lambda_index;

    for (ii=0; ii<db->pathNum; ii++)
    {
        if ((db->pathList[ii].wave[lambda-1] == 0) || (config->is_binary == 0))
        {
            path_num++;
        }
    }
    for (ii=0; ii<db->linkNum; ii++)
    {
        if ((db->linkList[ii].wave[lambda-1] == 0) ||  (config->is_binary == 0))
        {
            link_num++;
        }
    }
    if ((link_num == 0) || (path_num == 0))
    {
        for (ii=0; ii<db->pathNum; ii++)
        {
            config->config[ii] = 0;
        }
        return 0;
    }

    tmp_path_list = (RWA_PATH **)malloc(path_num * sizeof(RWA_PATH*));
    for (ii=0, jj=0; ii<db->pathNum; ii++)
    {
        if ((db->pathList[ii].wave[lambda-1] == 0) || (config->is_binary == 0))
        {
            tmp_path_list[jj++] = &db->pathList[ii];
        }
    }
    tmp_link_list = (RWA_LINK **)malloc(link_num * sizeof(RWA_LINK*));
    for (ii=0, jj=0; ii<db->linkNum; ii++)
    {
        if ((db->linkList[ii].wave[lambda-1] == 0) ||  (config->is_binary == 0))
        {
            tmp_link_list[jj++] = &db->linkList[ii];
        }
    }

    tmp_config = (double *)malloc(path_num * sizeof(double));

    *cost = d_pricing_problem_path(db, tmp_config, dual, link_num, path_num, tmp_link_list, tmp_path_list);

    for (ii=0; ii<path_num; ii++)
    {
        config->config[tmp_path_list[ii]->index-1] = tmp_config[ii];
    }

    free(tmp_link_list);
    free(tmp_path_list);
    free(tmp_config);

    return 1;
}

double d_pricing_problem_path(RWA_DB *db, double config[], const double dual[], const int link_num, const int path_num,
                              RWA_LINK *link_list[], RWA_PATH *path_list[])
{
    int ii,jj,kk;
    const int rows = link_num + db->reqGroupNum, columns = path_num;
    glp_prob *ppp;

    /* CAUTION!!! Every array data starts at element [1] instead of [0] - the element [0] belived to have its column or row's information. */
    int ia[1+LP_SIZE] = {0}, ja[1+LP_SIZE] = {0};
    double ar[1+LP_SIZE] = {0};
    double cost;

    ppp = glp_create_prob();
    glp_set_prob_name(ppp, "Pricing Problem_Path");
    glp_set_obj_dir(ppp, GLP_MAX);

    glp_add_rows(ppp, rows);
    for (ii = 1; ii < 1 + link_num; ii++)
    {
        glp_set_row_name(ppp, ii, "sum_p(delta(p,l)B[p])<=1");
        glp_set_row_bnds(ppp, ii, GLP_UP, 0.0, 1.0);
    }

    for (ii = 1 + link_num; ii < 1 + rows; ii++)//num_req_group
    {
        glp_set_row_name(ppp, ii, "sum_sd(B[p])<=D[sd]");
        glp_set_row_bnds(ppp, ii, GLP_UP, 0.0, db->reqList[ii - link_num - 1].req_num);
    }

    glp_add_cols(ppp, columns);
    for (ii = 1; ii < 1 + columns; ii++)//num_paths
    {
        //Beta[p]
        glp_set_col_name(ppp, ii, "B[p]");
        glp_set_col_kind(ppp, ii, GLP_BV);
        glp_set_obj_coef(ppp, ii, dual[path_list[ii-1]->req_index]);
    }

    for (ii = 1; ii < 1 + link_num; ii++)
    {
        for (jj = 1; jj < 1 + columns; jj++)//num_paths
        {
            ia[jj + columns * (ii - 1)] = ii;
            ja[jj + columns * (ii - 1)] = jj;
            for (kk = 1; kk < 1 + path_list[jj-1]->link_num; kk++)
            {
                if (path_list[jj-1]->link_list[kk-1] == link_list[ii-1]->index)
                {
                    ar[jj + columns * (ii - 1)] = 1; /* a[i,j] = delta(p,l) */
                    break;
                }
            }
        }
    }
    for (ii = 1 + link_num; ii < 1 + rows; ii++)//num_req_group
    {
        for (jj = 1; jj < 1 + columns; jj++)//num_paths
        {
            ia[jj + columns * (ii - 1)] = ii;
            ja[jj + columns * (ii - 1)] = jj;
            if (path_list[jj - 1]->req_index == ii - link_num)
            {
                ar[jj + columns * (ii - 1)] = 1; /* a[i,j] = delta(p, sd) */
            }
        }
    }

#ifdef PRINT_DETAIL_0
    printf("row bounds:\n");
    for (ii = 1; ii < 1 + link_num; ii++)
    {
        printf("sum_p(delta(p,l)B[%d])<=1\n", ii);
    }
    for (ii = 1 + link_num; ii < 1 + rows; ii++)//num_req_group
    {
        printf("sum_sd(B[p][%d])<=%d\n", ii - link_num - 1, db->reqList[ii - link_num - 1].req_num);
    }
    printf("var bounds:\n");
    for (ii = 1; ii < 1 + columns; ii++)//num_paths
    {
        printf("0<=B[%d]<=1\n", ii);
        printf("obj coef:req%d dual %f\n", db->reqList[path_list[ii-1]->req_index - 1].index, dual[path_list[ii-1]->req_index]);
    }
    printf("s.t. coef:\n");
    for (ii=1; ii<1+rows; ii++)
    {
        for (jj=1; jj<1+columns; jj++)
        {
            printf("%.1f\t", ar[jj + columns * (ii - 1)]);
        }
        printf("\n");
    }
#endif

    glp_load_matrix(ppp, rows * columns, ia, ja, ar);
    glp_simplex(ppp, NULL);
    glp_intopt(ppp, NULL);

    cost = glp_mip_obj_val(ppp) - dual[0];

    for (ii = 1; ii < 1 + columns; ii++)//num_paths
    {
        config[ii-1] = glp_mip_col_val(ppp, ii);
    }

    glp_delete_prob(ppp);

#ifdef PRINT_DETAIL_1
    printf("\nzPRICE = %g", cost);
    for (ii = 0; ii < path_num; ii++)
    {
        printf("\n path[%d] = %g", path_list[ii]->index, config[ii]);
    }
    printf("\n");
#endif

    return cost;
}

void test_dynamic_instance1()
{
    int link_list[MAX_SIZE] = {0};//parameter of init path; array of link index in a path.
    int path_list[MAX_SIZE] = {0};//parameter of init req. array of path index in a req.
    //int ii, jj;
    int W, linkNum, pathNum, reqGroupNum;
    RWA_LINK *linkList;
    RWA_PATH *pathList;
    RWA_REQ *reqList;
    RWA_DB *db;
    
    //allocate.
    db = (RWA_DB *)malloc(sizeof(RWA_DB));
    W = 10;
    linkNum = 6;
    pathNum = 6;
    reqGroupNum = 2;
    init_rwa_db(db, linkNum, pathNum, reqGroupNum, W);
    linkList = db->linkList;
    pathList = db->pathList;
    reqList = db->reqList;

    //init link data.
    init_link(&linkList[0], 1, 1, 4, W);
    init_link(&linkList[1], 2, 1, 2, W);
    init_link(&linkList[2], 3, 2, 3, W);
    init_link(&linkList[3], 4, 3, 4, W);
    init_link(&linkList[4], 5, 1, 3, W);
    init_link(&linkList[5], 6, 2, 4, W);

    //set wave status.
    set_link_wave(&linkList[2], 6, 1);
    set_link_wave(&linkList[2], 7, 1);
    set_link_wave(&linkList[2], 8, 1);
    set_link_wave(&linkList[2], 9, 1);
    set_link_wave(&linkList[2], 10, 1);
    //set_link_wave(&linkList[4], 4, 1);

    //init path list.
    link_list[0] = 5;
    init_path(&pathList[0], 1, 1, 3, 1, link_list, W);
    link_list[0] = 1; link_list[1] = 4;
    init_path(&pathList[1], 2, 1, 3, 2, link_list, W);
    link_list[0] = 2; link_list[1] = 3;
    init_path(&pathList[2], 3, 1, 3, 2, link_list, W);
    link_list[0] = 6;
    init_path(&pathList[3], 4, 2, 4, 1, link_list, W);
    link_list[0] = 2; link_list[1] = 1;
    init_path(&pathList[4], 5, 2, 4, 2, link_list, W);
    link_list[0] = 3; link_list[1] = 4;
    init_path(&pathList[5], 6, 2, 4, 2, link_list, W);
    //set_path_wave(pathList, pathNum, linkList, W);

    //init req data.
    path_list[0] = 1; path_list[1] = 2; path_list[2] = 3;
    init_req(&reqList[0], 1, 1, 3, 5, 3, path_list);
    path_list[0] = 4; path_list[1] = 5; path_list[2] = 6;
    init_req(&reqList[1], 2, 2, 4, 5, 3, path_list);

    set_path_req_index(reqList, reqGroupNum, pathList, pathNum);
    
    d_rwa(db);

    reset_path_result(pathList, pathNum, W);

    //destroy.
    del_db(db);
}

void test_dynamic_instance2()
{
    int link_list[MAX_SIZE] = {0};//parameter of init path; array of link index in a path.
    int path_list[MAX_SIZE] = {0};//parameter of init req. array of path index in a req.
    //int ii, jj;
    int W, linkNum, pathNum, reqGroupNum;
    RWA_LINK *linkList;
    RWA_PATH *pathList;
    RWA_REQ *reqList;
    RWA_DB *db;

    //allocate.
    db = (RWA_DB *)malloc(sizeof(RWA_DB));
    W = 10;
    linkNum = 6;
    pathNum = 3;
    reqGroupNum = 1;
    init_rwa_db(db, linkNum, pathNum, reqGroupNum, W);
    linkList = db->linkList;
    pathList = db->pathList;
    reqList = db->reqList;

    //init link data.
    init_link(&linkList[0], 1, 1, 4, W);
    init_link(&linkList[1], 2, 1, 2, W);
    init_link(&linkList[2], 3, 2, 3, W);
    init_link(&linkList[3], 4, 3, 4, W);
    init_link(&linkList[4], 5, 1, 3, W);
    init_link(&linkList[5], 6, 2, 4, W);

    //set wave status.
    //set_link_wave(&linkList[2], 2, 1);
    //set_link_wave(&linkList[4], 4, 1);

    //init path list.
    link_list[0] = 5;
    init_path(&pathList[0], 1, 1, 3, 1, link_list, W);
    link_list[0] = 1; link_list[1] = 4;
    init_path(&pathList[1], 2, 1, 3, 2, link_list, W);
    link_list[0] = 2; link_list[1] = 3;
    init_path(&pathList[2], 3, 1, 3, 2, link_list, W);
    //link_list[0] = 6;
    //init_path(&pathList[3], 4, 2, 4, 1, link_list, W);
    //link_list[0] = 2; link_list[1] = 1;
    //init_path(&pathList[4], 5, 2, 4, 2, link_list, W);
    //link_list[0] = 3; link_list[1] = 4;
    //init_path(&pathList[5], 6, 2, 4, 2, link_list, W);
    //set_path_wave(pathList, pathNum, linkList, W);

    //init req data.
    path_list[0] = 1; path_list[1] = 2; path_list[2] = 3;
    init_req(&reqList[0], 1, 1, 3, 10, 3, path_list);
    //path_list[0] = 4; path_list[1] = 5; path_list[2] = 6;
    //init_req(&reqList[1], 2, 2, 4, 10, 3, path_list);

    set_path_req_index(reqList, reqGroupNum, pathList, pathNum);

    d_rwa(db);

    reset_path_result(pathList, pathNum, W);

    //destroy.
    del_db(db);
}
