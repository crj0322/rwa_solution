/******************************************************************************************/
/*  rwa.c v_3.0 (MIS_based rwa solution) - Seyed Jafar Mosavian, Ruijia 10/8/2017         */
/******************************************************************************************/
#include "rwa.h"

void rwa(RWA_DB *db)
{
	//CAUTION!!! every array/matrix index starts from 0, except for the glpk arrays: ia[] ja[] ar[].
	int ii;
	int num_configurations = 0;
	double cost;
	double *dual;
	WAVE_CONFIG *config_list;
	int used_wave_num;
	int sum_ysd = 0;
	int sum_req = 0;
	int ret = 0;

	set_path_available_wavelength(db->path_list, db->path_num, db->link_list, db->W);

#ifdef PRINT_DETAIL_2
	print_req_list(db->req_list, db->req_qroup_num);
	print_path_list(db->path_list, db->path_num);
	print_link_list(db->link_list, db->link_num);
	print_path_wavelength(db->path_list, db->path_num, db->W);
	print_link_wavelenth(db->link_list, db->link_num, db->W);////-----------added
	//	print_path_result(db->pathList, db->path_num, db->W);
#endif

	//for (ii = 0; ii < db->req_qroup_num; ii++)//----------------????????????/---------------------------
	//{
	//    sum_req += db->req_list[ii].req_num;
	//}

	dual = (double *)malloc((db->req_qroup_num + 1) * sizeof(double));
	config_list = (WAVE_CONFIG *)malloc(db->W * sizeof(WAVE_CONFIG));
	used_wave_num = preset_config_list(db, config_list);
	//init dual values.
	dual[0] = 0;

	for (ii = 0; ii < db->req_qroup_num; ii++)
	{
		dual[ii + 1] = 1;
	}
	//dual[3] = 11;

	//solve.
	do
	{
		if (num_configurations == db->W)//------------------ *!!!*
		{
			break;
		}
		//generate wavelength configuration.
		config_list[num_configurations].config = (double *)malloc(db->path_num * sizeof(double));

		for (ii = 0; ii < db->path_num; ii++)
		{
			config_list[num_configurations].config[ii] = 0;
		}

		ret = column_generation(db, dual, &cost, &config_list[num_configurations]);

		//zeroing ac_sd[reqGroupNum].
		config_list[num_configurations].ac_sd = (int *)malloc(db->req_qroup_num * sizeof(int));

		for (ii = 0; ii < db->req_qroup_num; ii++)
		{
			config_list[num_configurations].ac_sd[ii] = 0;
		}

		if (ret == 0)
		{
			num_configurations++;// NOT available path for lambda = num_configurations - 1
			continue;
		}

		if (cost <= 0)//---------------??????????????????
		{
			break;
		}

		//calc ac_sd[reqGroupNum].
		for (ii = 0; ii < db->path_num; ii++)
		{
			if (config_list[num_configurations].config[ii])
			{
				config_list[num_configurations].ac_sd[db->path_list[ii].sub_req_index - 1]++;
			}
		}

#ifdef PRINT_DETAIL_0
		//print ac_sd.
		printf("ac_sd:\n");
		for (ii=0; ii<db->req_qroup_num; ii++)
		{
			for (jj=0; jj<=num_configurations; jj++)
			{
				printf("%d\t", config_list[jj].ac_sd[ii]);
			}
			printf("\n");
		}
#endif

		sum_ysd = master_problem(db, config_list, ++num_configurations, used_wave_num, dual);

	} while (cost > 0);// && (sum_ysd < sum_req));

	process_solutions(db, config_list, num_configurations, used_wave_num);

	//free configs.
	for (ii = 0; ii < db->W; ii++)
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
	//print_path_result(db->path_list, db->pathNum, db->W);////////////////////
	print_lightpath_list(db->lightpath_list, db->lightpath_num);
	print_path_wavelength(db->path_list, db->path_num, db->W);
	print_link_wavelenth(db->link_list, db->link_num, db->W);
	print_req_list(db->req_list, db->req_qroup_num);
#endif

	return;
}//end of rwa()

int preset_config_list(RWA_DB *db, WAVE_CONFIG config_list[])
{
	int ii, jj, kk;

	for (ii = 0; ii < db->W; ii++)
	{
		config_list[ii].is_binary = 0;
		config_list[ii].lambda_index = 0;
		config_list[ii].z_c = 0;
		config_list[ii].config = NULL;
		config_list[ii].ac_sd = NULL;
	}

	kk = 0;//index of config_list.

	for (ii = 0; ii < db->W; ii++)
	{
		//if any link had used this wavelength, add a binary z_c.
		for (jj = 0; jj < db->link_num; jj++)
		{
			if (db->link_list[jj].assigned_wavelength[ii])
			{
				config_list[kk].is_binary = 1;
				config_list[kk++].lambda_index = ii + 1;/////---------------------------------------   ??????? kk++  ?????
				break;
			}
		}
	}

	return kk;
}//end of preset_config_list()

void process_solutions(RWA_DB *db, const WAVE_CONFIG config_list[], const int num_configurations, const int used_wave_num)
{
	int ii, jj, kk;
	//	int lightpathIndex = 1;
	int start_wave_index;//start_wave_index for every configuration.
	int totalRequestNum;
	//	RWA_LINK *link;
	//	RWA_LIGHTPATH *lightpathList;
	int assignedWavelengthIndex[3];
	int *wave_status;//1-used by any link, 0-free for all links;
	int *wave_map;//index=index in config_list, value=db wave_index-1.

#ifdef PRINT_DETAIL_2
	//print solution.
	printf("solution:\n");
	printf("z_c\t");

	for (jj = 0; jj < num_configurations; jj++)
	{
		printf("%d\t", (int)config_list[jj].z_c);
	}

	printf("\n\n");

	for (ii = 0; ii < db->path_num; ii++)
	{
		printf("path[%d] ", ii + 1);

		for (jj = 0; jj < num_configurations; jj++)
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
	for (ii = 0; ii < db->W; ii++)
	{
		wave_status[ii] = 0;
	}

	for (ii = 0; ii < used_wave_num; ii++)
	{
		wave_map[ii] = config_list[ii].lambda_index;// - 1
		wave_status[wave_map[ii] - 1] = 1;
	}

	for (ii = 0, jj = used_wave_num; ii < db->W; ii++)
	{
		if (wave_status[ii] == 0)
		{
			wave_map[jj++] = ii + 1;// + 1
		}
	}

	printf("wave_map:");

	for (ii = 0; ii < db->W; ii++)
	{
		printf(" %d", wave_map[ii]);
	}

	printf("\n");

	if (jj != ii)
	{
		printf("calc wave_map error!\n");
	}

	for (ii = 0; ii < db->path_num; ii++)
	{
		for (jj = 0; jj < used_wave_num; jj++)
		{
			if ((config_list[jj].z_c) && (config_list[jj].config[ii]))
			{
				assignedWavelengthIndex[0] = config_list[jj].lambda_index ;// - 1
				//db->path_list[ii].assigned_wave_num++;/////////////////************************/////////////////
				db->lightpath_index++;
				init_lightpath(&db->lightpath_list[db->lightpath_index - 1], db->lightpath_index, 1, &db->path_list[ii], assignedWavelengthIndex);
				reset_lightpath_links(db, &db->lightpath_list[db->lightpath_index - 1]);
				//			db->lightpath_index++;
				//db->lightpathList[ii].assigned_wave_num++;
				//db->pathList[ii].assigned_wave[config_list[jj].lambda_index - 1] = 1;
				db->req_list[db->path_list[ii].sub_req_index - 1].req_num--;//-----------------------------------------------added-------------------
			}
		}

		start_wave_index = used_wave_num;//index of wave_map.
		totalRequestNum = db->req_list[db->path_list[ii].sub_req_index - 1].req_num;

		for (jj = used_wave_num; jj < num_configurations; jj++)
		{
			if ((config_list[jj].z_c) && (config_list[jj].config[ii]))
			{
				for (kk = 0; kk < (int)config_list[jj].z_c; kk++)
				{
					if (kk < totalRequestNum)
					{
					assignedWavelengthIndex[0] = wave_map[start_wave_index + kk];//******************************8
					//db->path_list[ii].assigned_wave_num++;
					//db->pathList[ii].assigned_wave[wave_map[start_wave_index + kk]] = 1;
					db->lightpath_index++;
					init_lightpath(&db->lightpath_list[db->lightpath_index - 1], db->lightpath_index, 1, &db->path_list[ii], assignedWavelengthIndex);
					reset_lightpath_links(db, &db->lightpath_list[db->lightpath_index - 1]);
					//	lightpathIndex++;
					//					db->reqList[db->pathList[ii].req_index - 1].req_num--;//-----------------------------------------------added-------------------

					db->req_list[db->path_list[ii].sub_req_index - 1].req_num--;//-----------------------------------------------added-------------------
					}
				}
			}

			if (config_list[jj].z_c)
			{
				start_wave_index += (int)(config_list[jj].z_c);
			}
		}

		//set wave status to linkList.
		//	for (jj = 0; jj < db->pathList[ii].link_num; jj++)
		//	{
		//		link = &db->linkList[db->pathList[ii].link_list[jj] - 1];

		//		if (db->pathList[ii].assigned_wave_num > 0)
		//		{
		//			for (kk = 0; kk < db->W; kk++)
		//			{
		//				if (link->wave[kk] && db->pathList[ii].assigned_wave[kk])
		//				{
		//					printf("path%d link%d wave%d conflict!!!\n", ii + 1, db->pathList[ii].link_list[jj], kk + 1);
		//					//process conflict.
		//					//to do...
		//				}

		//				link->wave[kk] |= db->pathList[ii].assigned_wave[kk];
		//			}
		//		}
		//	}
	}

	free(wave_status);
	free(wave_map);
	//set_path_wave(db->pathList, db->pathNum, db->linkList, db->W);
	return;
}//end of process_solutions()

int master_problem(RWA_DB *db, WAVE_CONFIG config_list[], const int num_configurations, const int W_used, double dual[])
{
	int ii, jj;
	glp_prob *mp;
	const int rows = 1 + db->req_qroup_num, columns = num_configurations + db->req_qroup_num;
	double obj;
	int sum_ysd = 0;
	/* CAUTION!!! Every array data starts at element [1] instead of [0] - the element [0] belived to have its column or row's information. */
	int ia[1 + LP_SIZE] = {0}, ja[1 + LP_SIZE] = {0};
	double ar[1 + LP_SIZE] = {0};

	mp = glp_create_prob();

	glp_set_prob_name(mp, "(Restricted) Master Problem");
	glp_set_obj_dir(mp, GLP_MAX);
	glp_add_rows(mp, rows);

	glp_set_row_name(mp, 1, "sum_c(z_c)<=W");
	glp_set_row_bnds(mp, 1, GLP_UP, 0.0, db->W - W_used);

	for (ii = 2; ii < 2 + db->req_qroup_num; ii++)
	{
		glp_set_row_name(mp, ii, "y[sd]-sum_c(ac[sd]*z_c)<=0");
		glp_set_row_bnds(mp, ii, GLP_UP, 0.0, 0.0);
	}

	glp_add_cols(mp, columns);

	for (ii = 1; ii < 1 + db->req_qroup_num; ii++)
	{
		glp_set_col_name(mp, ii, "y_sd");
		glp_set_col_bnds(mp, ii, GLP_DB, 0.0, db->req_list[ii - 1].req_num);
		glp_set_col_kind(mp, ii, GLP_IV);
		glp_set_obj_coef(mp, ii, 1.0);
	}

	for (ii = 1 + db->req_qroup_num; ii < 1 + columns; ii++)//num_configurations
	{
		glp_set_col_name(mp, ii, "z_c");

		if (config_list[ii - db->req_qroup_num - 1].is_binary)
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

	for (jj = 1; jj < 1 + db->req_qroup_num; jj++)
	{
		ia[jj + columns * (ii - 1)] = ii;
		ja[jj + columns * (ii - 1)] = jj;
		ar[jj + columns * (ii - 1)] = 0; // a[i,j] = 0 for y_sd.
	}

	for (jj = 1 + db->req_qroup_num; jj < 1 + columns; jj++)//num_configurations
	{
		ia[jj + columns * (ii - 1)] = ii;
		ja[jj + columns * (ii - 1)] = jj;
		ar[jj + columns * (ii - 1)] = 1; // a[i,j] = 1 for z_c.

		if (config_list[jj - db->req_qroup_num - 1].is_binary)//---------------------------added
		{
			ar[jj + columns * (ii - 1)] = 0; // a[i,j] = 0 for z_c used.
		}
	}

	for (ii = 2; ii < 2 + db->req_qroup_num; ii++)
	{
		for (jj = 1; jj < 1 + db->req_qroup_num; jj++)
		{
			ia[jj + columns * (ii - 1)] = ii;
			ja[jj + columns * (ii - 1)] = jj;

			if (jj == ii - 1)
			{
				ar[jj + columns * (ii - 1)] = 1; // a[i,j] = 1 for y_sd[jj].
			}
		}

		for (jj = 1 + db->req_qroup_num; jj < 1 + columns; jj++)//num_configurations
		{
			ia[jj + columns * (ii - 1)] = ii;
			ja[jj + columns * (ii - 1)] = jj;
			ar[jj + columns * (ii - 1)] = -(double)config_list[jj - db->req_qroup_num - 1].ac_sd[ii - 2]; // a[i,j] = -ac_sd[i,j].
		}
	}

#ifdef PRINT_DETAIL_0
	printf("row bounds:\n");
	printf("sum_c(z_c)<=W\n");
	for (ii = 2; ii < 2 + db->req_qroup_num; ii++)
	{
		printf("y_sd[%d]-sum_c(ac[%d]*z_c)<=0\n", ii - 1, ii - 1);
	}
	printf("var bounds:\n");
	for (ii = 1; ii < 1 + db->req_qroup_num; ii++)
	{
		printf("0<=y_sd[%d]<=%d\n", ii, db->req_list[ii - 1].req_num);
	}
	for (ii = 1 + db->req_qroup_num; ii < 1 + columns; ii++)//num_configurations
	{
		printf("z_c[%d]>=0\n", ii - db->req_qroup_num);
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

	for (ii = 1; ii < 1 + db->req_qroup_num; ii++)
	{
		sum_ysd += (int)glp_mip_col_val(mp, ii);
	}

	for (ii = 1 + db->req_qroup_num; ii < 1 + columns; ii++)//num_configurations
	{
		config_list[ii - db->req_qroup_num - 1].z_c = glp_mip_col_val(mp, ii);
	}

	//get dual values.
	dual[0] = glp_get_row_dual(mp, 1);

	for (ii = 2; ii < 2 + db->req_qroup_num; ii++)
	{
		dual[ii - 1] = glp_get_row_dual(mp, ii);
	}

#ifdef PRINT_DETAIL_1
	printf("\nMASTER obj = %g", obj);
	for (ii = 1; ii < 1 + db->req_qroup_num; ii++)
	{
		printf("\n y[%d] = %g", ii, glp_mip_col_val(mp, ii));
	}
	printf("\n");
	for (ii = 1 + db->req_qroup_num; ii < 1 + columns; ii++)//num_configurations
	{
		printf("\n z_c[%d] = %g", ii - db->req_qroup_num, config_list[ii - db->req_qroup_num - 1].z_c);
	}
	printf("\n");
#endif

	return sum_ysd;
}//end of master_problem()

int column_generation(RWA_DB *db, const double dual[], double *cost, WAVE_CONFIG *config)
{
	int link_num = 0, path_num = 0, ii, jj;
	RWA_LINK **tmp_link_list;
	RWA_PATH **tmp_path_list;
	double *tmp_config, cst;
	int lambda = config->lambda_index;

	for (ii = 0; ii < db->path_num; ii++)
	{
		if ((db->path_list[ii].available_wavelength[lambda - 1] == 0) || (config->is_binary == 0))
		{
			path_num++;
		}
	}

	for (ii = 0; ii < db->link_num; ii++)
	{
		if ((db->link_list[ii].assigned_wavelength[lambda - 1] == 0) ||  (config->is_binary == 0))
		{
			link_num++;
		}
	}

	if ((link_num == 0) || (path_num == 0))
	{
		for (ii = 0; ii < db->path_num; ii++)
		{
			config->config[ii] = 0;
		}
		*cost = 800; // this means the current lamda is skipped from present solution process
		return 0;
	}

	tmp_path_list = (RWA_PATH **)malloc(path_num * sizeof(RWA_PATH*));

	for (ii = 0, jj = 0; ii < db->path_num; ii++)
	{
		if ((db->path_list[ii].available_wavelength[lambda - 1] == 0) || (config->is_binary == 0))
		{
			tmp_path_list[jj++] = &db->path_list[ii];
		}
	}

	tmp_link_list = (RWA_LINK **)malloc(link_num * sizeof(RWA_LINK*));

	for (ii = 0, jj = 0; ii < db->link_num; ii++)
	{
		if ((db->link_list[ii].assigned_wavelength[lambda - 1] == 0) ||  (config->is_binary == 0))
		{
			tmp_link_list[jj++] = &db->link_list[ii];
		}
	}

	tmp_config = (double *)malloc(path_num * sizeof(double));

	cst = pricing_problem_path(db, tmp_config, dual, link_num, path_num, tmp_link_list, tmp_path_list);

	for (ii = 0; ii < path_num; ii++)
	{
		config->config[tmp_path_list[ii]->index - 1] = tmp_config[ii];
	}
	///--------- added

	*cost = cst;

	if (cst <= 0)
	{
		for (ii = 1; ii < 1 + db->req_qroup_num; ii++)
		{
			if (dual[ii] > 0)
			{
				*cost = 900; // this is when present lamda (configuration) is not produced any solution!
				break;
			}
		}
	}

	///----------

	free(tmp_link_list);
	free(tmp_path_list);
	free(tmp_config);
	return 1;
}//end of column_generation()

double pricing_problem_path(RWA_DB *db, double config[], const double dual[], const int link_num, const int path_num,
							RWA_LINK *link_list[], RWA_PATH *path_list[])
{
	int ii, jj, kk;
	const int rows = link_num + db->req_qroup_num, columns = path_num;
	glp_prob *ppp;

	/* CAUTION!!! Every array data starts at element [1] instead of [0] - the element [0] belived to have its column or row's information. */
	int ia[1 + LP_SIZE] = {0}, ja[1 + LP_SIZE] = {0};
	double ar[1 + LP_SIZE] = {0};
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
		glp_set_row_bnds(ppp, ii, GLP_UP, 0.0, db->req_list[ii - link_num - 1].req_num);
	}

	glp_add_cols(ppp, columns);

	for (ii = 1; ii < 1 + columns; ii++)//num_paths
	{
		//Beta[p]
		glp_set_col_name(ppp, ii, "B[p]");
		glp_set_col_kind(ppp, ii, GLP_BV);
		glp_set_obj_coef(ppp, ii, dual[path_list[ii - 1]->sub_req_index]);
	}

	for (ii = 1; ii < 1 + link_num; ii++)
	{
		for (jj = 1; jj < 1 + columns; jj++)//num_paths
		{
			ia[jj + columns * (ii - 1)] = ii;
			ja[jj + columns * (ii - 1)] = jj;
			for (kk = 1; kk < 1 + path_list[jj - 1]->link_num; kk++)
			{
				if (path_list[jj - 1]->link_list[kk - 1] == link_list[ii - 1]->index)
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
			if (path_list[jj - 1]->sub_req_index == ii - link_num)
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
		printf("sum_sd(B[p][%d])<=%d\n", ii - link_num - 1, db->req_list[ii - link_num - 1].req_num);
	}
	printf("var bounds:\n");
	for (ii = 1; ii < 1 + columns; ii++)//num_paths
	{
		printf("0<=B[%d]<=1\n", ii);
		printf("obj coef:req%d dual %f\n", db->req_list[path_list[ii-1]->req_index - 1].index, dual[path_list[ii-1]->req_index]);
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
}//end of pricing_problem_path()


void test_instance_3()
{
	int link_list[MAX_SIZE] = {0};//parameter of init path; array of link index in a path.
	int path_list[MAX_SIZE] = {0};//parameter of init req. array of path index in a req.
	int ww;
	int W, linkNum, pathNum, reqGroupNum;
	int lightpathNum = 150;//////////////////////////////////// **********************/////////////////////////
	RWA_LINK *linkList;
	RWA_PATH *pathList;
	//	RWA_LIGHTPATH *lightpathList;
	RWA_REQ *reqList;
	RWA_DB *db;

	//allocate.
	db = (RWA_DB *)malloc(sizeof(RWA_DB));
	W = 80;
	linkNum = 10;
	pathNum = 5;
	reqGroupNum = 4;
	//	init_rwa_db(db, linkNum, pathNum, reqGroupNum, W);
	init_rwa_db(db, linkNum, pathNum, lightpathNum, reqGroupNum, W);
	linkList = db->link_list;
	pathList = db->path_list;
	//	lightpathList = db->lightpathList;
	reqList = db->req_list;

	//init link data.
	init_link(&linkList[0], 1, 1, 2, W);
	init_link(&linkList[1], 2, 1, 3, W);
	init_link(&linkList[2], 3, 4, 2, W);
	init_link(&linkList[3], 4, 2, 5, W);
	init_link(&linkList[4], 5, 3, 4, W);
	init_link(&linkList[5], 6, 5, 4, W);
	init_link(&linkList[6], 7, 3, 6, W);
	init_link(&linkList[7], 8, 4, 7, W);
	init_link(&linkList[8], 9, 5, 7, W);
	init_link(&linkList[9], 10, 6, 7, W);

	//set wavelength status.
	for (ww = 1; ww <= 80; ww++)
		set_link_wavelength(&linkList[7], ww, 1);
	for (ww = 1; ww <= 30; ww++)
		set_link_wavelength(&linkList[3], ww, 1);
	for (ww = 41; ww <= 75; ww++)
		set_link_wavelength(&linkList[3], ww, 1);
	for (ww = 1; ww <= 10; ww++)
		set_link_wavelength(&linkList[2], ww, 1);
	for (ww = 21; ww <= 30; ww++)
		set_link_wavelength(&linkList[2], ww, 1);
	for (ww = 11; ww <= 18; ww++)
		set_link_wavelength(&linkList[4], ww, 1);
	for (ww = 31; ww <= 78; ww++)
		set_link_wavelength(&linkList[4], ww, 1);
	for (ww = 1; ww <= 30; ww++)
		set_link_wavelength(&linkList[6], ww, 1);
	for (ww = 11; ww <= 20; ww++)
		set_link_wavelength(&linkList[8], ww, 1);
	for (ww = 31; ww <= 40; ww++)
		set_link_wavelength(&linkList[8], ww, 1);


	//init path list.
	link_list[0] = 1; link_list[1] = 4;
	init_path(&pathList[0], 1, 1, 5, 2, link_list, W);
	link_list[0] = 6; link_list[1] = 9;
	init_path(&pathList[1], 2, 4, 7, 2, link_list, W);
	link_list[0] = 3; link_list[1] = 4; link_list[2] = 9;
	init_path(&pathList[2], 3, 4, 7, 3, link_list, W);
	link_list[0] = 4; link_list[1] = 9;
	init_path(&pathList[3], 4, 2, 7, 2, link_list, W);
	link_list[0] = 10; link_list[1] = 9; link_list[2] = 4;
	init_path(&pathList[4], 5, 6, 2, 3, link_list, W);

	//	set_path_available_wavelength(pathList, pathNum, linkList, W);

	//init req data.
	path_list[0] = 1;
	init_req(&reqList[0], 1, 1, 5, 20, 1, path_list);
	path_list[0] = 2; path_list[1] = 3;
	init_req(&reqList[1], 2, 4, 7, 60, 2, path_list);
	path_list[0] = 4;
	init_req(&reqList[2], 3, 2, 7, 10, 1, path_list);
	path_list[0] = 5;
	init_req(&reqList[3], 4, 6, 2, 10, 1, path_list);

	set_path_req_index(reqList, reqGroupNum, pathList, pathNum);

	rwa(db);

	printf("\n\n\n\nfirst group is done!\n\n\n");

	reset_rwa_db(db, pathNum, reqGroupNum);

	pathList = db->path_list;
	reqList = db->req_list;




	//init path list.
	link_list[0] = 1; link_list[1] = 4;
	init_path(&pathList[0], 1, 1, 5, 2, link_list, W);
	link_list[0] = 6; link_list[1] = 9;
	init_path(&pathList[1], 2, 4, 7, 2, link_list, W);
	link_list[0] = 3; link_list[1] = 4; link_list[2] = 9;
	init_path(&pathList[2], 3, 4, 7, 3, link_list, W);
	link_list[0] = 4; link_list[1] = 9;
	init_path(&pathList[3], 4, 2, 7, 2, link_list, W);
	link_list[0] = 10; link_list[1] = 9; link_list[2] = 4;
	init_path(&pathList[4], 5, 6, 2, 3, link_list, W);

	//	set_path_available_wavelength(pathList, pathNum, linkList, W);

	//init req data.
	path_list[0] = 1;
	init_req(&reqList[0], 1, 1, 5, 10, 1, path_list);
	path_list[0] = 2; path_list[1] = 3;
	init_req(&reqList[1], 2, 4, 7, 10, 2, path_list);
	path_list[0] = 4;
	init_req(&reqList[2], 3, 2, 7, 10, 1, path_list);
	path_list[0] = 5;
	init_req(&reqList[3], 4, 6, 2, 10, 1, path_list);

	set_path_req_index(reqList, reqGroupNum, pathList, pathNum);

	rwa(db);

	printf("\n\n\n\n reqGroupNum = %d\n\n\n",db->req_qroup_num);


	del_db(db);

}//end of test_instance_3()
