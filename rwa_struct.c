/******************************************************************************************/
/*  rwa_struct.c v_3.0 (MIS_based rwa solution) - Seyed Jafar Mosavian, Ruijia 10/8/2017  */
/******************************************************************************************/
#include "head.h"
#include "rwa_struct.h"

void init_link(RWA_LINK *link, int index, int src, int des, int W)
{
    int ii;
    link->index = index;
    link->src = src;
    link->des = des;
    link->assigned_wavelength = (int *)malloc(W * sizeof(int));

    for (ii = 0; ii < W; ii++)
    {
        link->assigned_wavelength[ii] = 0;//init as all free.
    }

	link->lightpath_list = (int *)malloc(W * sizeof(int));

	for (ii = 0; ii < W; ii++)
    {
        link->lightpath_list[ii] = 0;//init as all free.
    }
}

void print_link_list(RWA_LINK *linkList, int linkNum)
{
    int ii;
    RWA_LINK *link;

    for (ii = 0; ii < linkNum; ii++)
    {
        link = &linkList[ii];
        printf("link%d(S-D:%d-%d)\n", link->index, link->src, link->des);
    }

    printf("\n");
}

void del_link_list(RWA_LINK *linkList, int linkNum)
{
    int ii;

    for (ii = 0; ii < linkNum; ii++)
    {
        if (linkList[ii].assigned_wavelength != NULL)
        {
            free(linkList[ii].assigned_wavelength);
        }

		if (linkList[ii].lightpath_list != NULL)
        {
            free(linkList[ii].lightpath_list);//before freeing memory, we need to record lightpaths of the link to our new requests
        }
    }
    
    free(linkList);
    linkList = NULL;
    linkNum = 0;
}

void set_link_wavelength(RWA_LINK *link, int wave_index, int wave_status)
{
    link->assigned_wavelength[wave_index - 1] = wave_status;
}

void print_link_wavelenth(RWA_LINK *linkList, int linkNum, int W)
{
    int ii, jj;
    printf("link wavelength info:\n");

    for (ii = 0; ii < linkNum; ii++)
    {
        printf("link[%d]:", linkList[ii].index);

        for (jj = 0; jj < W; jj++)
        {
            printf(" %d", linkList[ii].assigned_wavelength[jj]);
        }

        printf("\n");
    }

    printf("\n");
}

void init_path(RWA_PATH *path, int index, int src, int des, int link_num, int link_list[], int W)
{
    int ii;
    path->index = index;
    path->src = src;
    path->des = des;
    path->link_num = link_num;
    path->link_list = (int*)malloc(link_num * sizeof(int));

    for (ii = 0; ii < link_num; ii++)
    {
        path->link_list[ii] = link_list[ii];
    }

	path->sub_req_index = 0;
	path->available_wavelength = (int *)malloc(W * sizeof(int));

	for (ii = 0; ii < W; ii++)
    {
        path->available_wavelength[ii] = 0;//init as all free.
    }
}

void print_path_list(RWA_PATH *pathList, int pathNum)
{
    int ii, jj;
    RWA_PATH *path;

    for (jj = 0; jj < pathNum; jj++)
    {
        path = &pathList[jj];
        printf("path%d(S-D:%d-%d)(link", path->index, path->src, path->des);

        for (ii = 0; ii < path->link_num; ii++)
        {
            printf("-%d", path->link_list[ii]);
        }

        printf(")\n");
    }

    printf("\n");
}

void del_path_list(RWA_PATH *pathList, int pathNum)
{
    int ii;

    for (ii = 0; ii < pathNum; ii++)
    {
        if (pathList[ii].link_list != NULL)
        {
            free(pathList[ii].link_list);
            pathList[ii].link_list = NULL;
        }
	}

	for (ii = 0; ii < pathNum; ii++)
	{
		if (pathList[ii].available_wavelength != NULL)
		{
			free(pathList[ii].available_wavelength);
			pathList[ii].available_wavelength = NULL;
		}
	}

	free(pathList);
    pathList = NULL;
    pathNum = 0;
}

void set_path_available_wavelength(RWA_PATH pathList[], int pathNum, RWA_LINK linkList[], int W)
{
	int ii, jj, kk;

	for (ii = 0; ii < pathNum; ii++)
	{
		for (jj = 0; jj < W; jj++)
		{
			pathList[ii].available_wavelength[jj] = 0;

			for (kk = 0; kk < pathList[ii].link_num; kk++)
			{
				pathList[ii].available_wavelength[jj] |= linkList[pathList[ii].link_list[kk] - 1].assigned_wavelength[jj];
			}
		}
	}
}

void print_path_wavelength(RWA_PATH *pathList, int pathNum, int W)
{
    int ii, jj;
    printf("path wavelength info:\n");

    for (ii = 0; ii < pathNum; ii++)
    {
        printf("path[%d]:", pathList[ii].index);

        for (jj = 0; jj < W; jj++)
        {
            printf(" %d", pathList[ii].available_wavelength[jj]);
        }

        printf("\n");
    }

    printf("\n");
}

void init_lightpath(RWA_LIGHTPATH *lightpath, int index, int subPathNum, RWA_PATH subPathList[], int assignedWavelengthIndex[])
{
    int ii;
//	db->lightpath_index++;
    lightpath->index = index;
    lightpath->sub_path_num = subPathNum;
	lightpath->sub_path_list = (RWA_PATH*)malloc(subPathNum * sizeof(RWA_PATH));

    for (ii = 0; ii < subPathNum; ii++)
    {
        lightpath->sub_path_list[ii] = subPathList[ii];
    }

    lightpath->assigned_wavelength_index = (int*)malloc(subPathNum * sizeof(int));

    for (ii = 0; ii < subPathNum; ii++)
    {
        lightpath->assigned_wavelength_index[ii] = assignedWavelengthIndex[ii];
    }

	//lightpath->req_index = reqIndex;
}

void print_lightpath_list(RWA_LIGHTPATH *lightpathList, int lightpathNum)
{
	int ii, jj;
	printf("lightpath info:\n");

	for (ii = 0; ii < lightpathNum; ii++)
	{
		printf("lightpath[%d]:", lightpathList[ii].index);

		for (jj = 0; jj < lightpathList[ii].sub_path_num; jj++)
		{
			printf(" S-D %d-%d path_index: %d", lightpathList[ii].sub_path_list[jj].src, lightpathList[ii].sub_path_list[jj].des, lightpathList[ii].sub_path_list[jj].index);
			printf(" assigned_lambda %d", lightpathList[ii].assigned_wavelength_index[jj]);
		}

		printf("\n");
	}

	printf("\n");
}

void del_lightpath_list(RWA_LIGHTPATH *lightpathList, int lightpathIndex)
{
	int ii;

	for (ii = 0; ii < lightpathIndex; ii++)
	{
		/*for (jj = 0; jj < lightpathList[ii].sub_path_num; jj++)
		{
			if (lightpathList[ii].sub_path_list[jj] != NULL)
			{
				free(lightpathList[ii].sub_path_list[jj]);
				lightpathList[ii].sub_path_list[jj] = NULL;
			}
		}*/
		free(lightpathList[ii].sub_path_list);
		free(lightpathList[ii].assigned_wavelength_index);
		/*for (kk = 0; kk < lightpathList[ii].sub_path_num; kk++)
		{
			if (lightpathList[ii].assigned_wavelength[kk] != NULL)
			{
				free(lightpathList[ii].assigned_wavelength[kk]);
				lightpathList[ii].assigned_wavelength[kk] = NULL;
			}
		}*/
	}

	free(lightpathList);
	lightpathList = NULL;
	lightpathIndex = 0;
}

void reset_lightpath_links(RWA_DB *db, RWA_LIGHTPATH *lightpath)
{
	int ii, jj;
	RWA_LINK *link;

	for (ii = 0; ii < lightpath->sub_path_num; ii++)
	{
		if (lightpath->assigned_wavelength_index[ii] > 0)
		{
			for (jj = 0; jj < lightpath->sub_path_list[ii].link_num; jj++)
			{
				link = &db->link_list[lightpath->sub_path_list[ii].link_list[jj] - 1];
				link->assigned_wavelength[lightpath->assigned_wavelength_index[ii] - 1] = 1;//assigned_wavelength[lightpath->assigned_wavelength_index[ii] - 1]
			}
		}
	}
}

void init_req(RWA_REQ *req, int index, int src, int des, int reqNum, int pathNum, int pathList[])
{
    int ii;
    req->index = index;
    req->src = src;
    req->des = des;
    req->path_num = pathNum;
    req->req_num = reqNum;
    req->path_list = (int *)malloc(pathNum * sizeof(int));

    for (ii = 0; ii < pathNum; ii++)
    {
        req->path_list[ii] = pathList[ii];
    }
}

void del_req_list(RWA_REQ *reqList, int reqGroupNum)
{
    int ii;

    for (ii = 0; ii < reqGroupNum; ii++)
    {
        if (reqList[ii].path_list != NULL)
        {
            free(reqList[ii].path_list);
        }
    }

    free(reqList);
    reqList = NULL;
    reqGroupNum = 0;
}

void print_req_list(RWA_REQ *reqList, int reqGroupNum)
{
    int ii, jj;
    RWA_REQ *req;

    for (jj = 0; jj < reqGroupNum; jj++)
    {
        req = &reqList[jj];
        printf("req%d(SD:%d-%d)(req_num:%d)(path", req->index, req->src, req->des, req->req_num);

        for (ii = 0; ii < req->path_num; ii++)
        {
            printf(" %d", req->path_list[ii]);
        }

        printf(")\n");
    }

    printf("\n");
}

void set_path_req_index(RWA_REQ *reqList, int reqGroupNum, RWA_PATH *pathList, int pathNum)
{
    int ii, jj;
    RWA_PATH *path;
    RWA_REQ *req;
    
    for (ii = 0; ii < pathNum; ii++)
    {
        path = &pathList[ii];

        for (jj = 0; jj < reqGroupNum; jj++)
        {
            req = &reqList[jj];

            if ((path->src == req->src) && (path->des == req->des))
            {
                path->sub_req_index = req->index;
            }
        }
    }
}

void init_rwa_db(RWA_DB *db, int linkNum, int pathNum, int lightpathNum, int reqGroupNum, int W)
{
    db->W = W;
	db->link_num = linkNum;
    db->path_num = pathNum;
    db->req_qroup_num = reqGroupNum;
	db->lightpath_index = 0;
	db->lightpath_num = lightpathNum;
    db->req_list = (RWA_REQ*)malloc(reqGroupNum * sizeof(RWA_REQ));
    db->path_list = (RWA_PATH*)malloc(pathNum * sizeof(RWA_PATH));
    db->link_list = (RWA_LINK*)malloc(linkNum * sizeof(RWA_LINK));
	db->lightpath_list = (RWA_LIGHTPATH*)malloc(lightpathNum * sizeof(RWA_LIGHTPATH));//  needs more work!  ***lightpathNum is temporary! ****
}

void reset_rwa_db(RWA_DB *db, int pathNum, int reqGroupNum)
{
	del_path_list(db->path_list, db->path_num);
	del_req_list(db->req_list, db->req_qroup_num);
	//    db->W = W;
	//	db->linkNum = linkNum;
	db->path_num = pathNum;
    db->req_qroup_num = reqGroupNum;
//	db->lightpath_index = 0;
//	db->lightpathNum = lightpathNum;
    db->req_list = (RWA_REQ*)malloc(reqGroupNum * sizeof(RWA_REQ));
    db->path_list = (RWA_PATH*)malloc(pathNum * sizeof(RWA_PATH));
//    db->link_list = (RWA_LINK*)malloc(linkNum * sizeof(RWA_LINK));
//	db->lightpathList = (RWA_LIGHTPATH*)malloc(lightpathNum * sizeof(RWA_LIGHTPATH));//  needs more work!  ***lightpathNum is temporary! ****
}

void del_db(RWA_DB *db)
{
    if (db != NULL)
    {
        del_link_list(db->link_list, db->link_num);
        del_path_list(db->path_list, db->path_num);
        del_req_list(db->req_list, db->req_qroup_num);
		del_lightpath_list(db->lightpath_list, db->lightpath_index);
        free(db);
    }
}