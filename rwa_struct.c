/******************************************************************************************/
/*  rwa_struct.c v_1.0 (MIS_based rwa solution) - Seyed Jafar Mosavian, Ruijia 9/11/2017  */
/******************************************************************************************/
#include "head.h"
#include "rwa_struct.h"

void init_link(RWA_LINK *link, int index, int src, int des, int W)
{
    int ii;
    link->index = index;
    link->src = src;
    link->des = des;
    link->wave = (int *)malloc(W * sizeof(int));
    for (ii=0; ii<W; ii++)
    {
        link->wave[ii] = 0;//init as all free.
    }
}

void print_link_list(RWA_LINK *linkList, int linkNum)
{
    int ii;
    RWA_LINK *link;
    for (ii=0; ii<linkNum; ii++)
    {
        link = &linkList[ii];
        printf("link%d(SD:%d-%d)\n", link->index, link->src, link->des);
    }
    printf("\n");
}

void del_link_list(RWA_LINK *linkList, int linkNum)
{
    int ii;

    for (ii=0; ii<linkNum; ii++)
    {
        if (linkList[ii].wave != NULL)
        {
            free(linkList[ii].wave);
        }
    }
    
    free(linkList);
    linkList = NULL;
    linkNum = 0;
}

void print_link_wave(RWA_LINK *linkList, int linkNum, int W)
{
    int ii, jj;

    printf("link wave info:\n");
    for (ii = 0; ii < linkNum; ii++)
    {
        printf("link[%d]:", linkList[ii].index);
        for (jj = 0; jj < W; jj++)
        {
            printf(" %d", linkList[ii].wave[jj]);
        }
        printf("\n");
    }
    printf("\n");
}

void set_link_wave(RWA_LINK *link, int wave_index, int wave_status)
{
    link->wave[wave_index - 1] = wave_status;
}

void init_path(RWA_PATH *path, int index, int src, int des, int link_num, int link_list[], int W)
{
    int ii;
    path->index = index;
    path->src = src;
    path->des = des;
    path->link_num = link_num;
    path->req_index = 0;
    path->link_list = (int*)malloc(link_num * sizeof(int));
    for (ii=0; ii<link_num; ii++)
    {
        path->link_list[ii] = link_list[ii];
    }
    path->wave = (int *)malloc(W * sizeof(int));
    for (ii = 0; ii < W; ii++)
    {
        path->wave[ii] = 0;
    }

    path->assigned_wave_num = 0;
    path->assigned_wave = (int *)malloc(W * sizeof(int));
    for (ii = 0; ii < W; ii++)
    {
        path->assigned_wave[ii] = 0;
    }
}

void print_path_list(RWA_PATH *pathList, int pathNum)
{
    int ii, jj;
    RWA_PATH *path;
    for (jj = 0; jj < pathNum; jj++)
    {
        path = &pathList[jj];
        printf("path%d(SD:%d-%d)(link", path->index, path->src, path->des);
        for (ii=0; ii<path->link_num; ii++)
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
        if (pathList[ii].wave != NULL)
        {
            free(pathList[ii].wave);
            pathList[ii].wave = NULL;
        }
        if (pathList[ii].assigned_wave != NULL)
        {
            free(pathList[ii].assigned_wave);
            pathList[ii].assigned_wave = NULL;
        }
    }

    free(pathList);
    pathList = NULL;
    pathNum = 0;
}

void set_path_wave(RWA_PATH pathList[], int pathNum, RWA_LINK linkList[], int W)
{
	int ii, jj, kk;
	for (ii = 0; ii < pathNum; ii++)
	{
		for (jj = 0; jj < W; jj++)
		{
			//pathList[ii].wave[jj] = 0;
			for (kk = 0; kk < pathList[ii].link_num; kk++)
			{
				pathList[ii].wave[jj] |= linkList[pathList[ii].link_list[kk] - 1].wave[jj];
			}
		}
	}
}

void reset_path_result(RWA_PATH pathList[], int pathNum, int W)
{
    int ii, jj;
    for (ii = 0; ii < pathNum; ii++)
    {
        pathList[ii].assigned_wave_num = 0;
        if (pathList[ii].assigned_wave)
        {
            for (jj = 0; jj < W; jj++)
            {
                pathList[ii].assigned_wave[jj] = 0;
            }
        }
    }
}

void print_path_wave(RWA_PATH *pathList, int pathNum, int W)
{
    int ii, jj;

    printf("path wave info:\n");
    for (ii = 0; ii < pathNum; ii++)
    {
        printf("path[%d]:", pathList[ii].index);
        for (jj = 0; jj < W; jj++)
        {
            printf(" %d", pathList[ii].wave[jj]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_path_result(RWA_PATH pathList[], int pathNum, int W)
{
    int ii, jj;
    printf("path result:\n");
    for (ii = 0; ii < pathNum; ii++)
    {
        if ((pathList[ii].assigned_wave_num > 0) && (pathList[ii].assigned_wave != NULL))
        {
            printf("assigned wave of path[%d]: ", ii + 1);
            for (jj = 0; jj < W; jj++)
            {
                printf("%d ", pathList[ii].assigned_wave[jj]);
            }
            printf("\n");
        }
    }
    printf("\n");
}

void init_req(RWA_REQ *req, int index, int src, int des, int req_num, int path_num, int path_list[])
{
    int ii;
    req->index = index;
    req->src = src;
    req->des = des;
    req->path_num = path_num;
    req->req_num = req_num;
    req->path_list = (int *)malloc(path_num * sizeof(int));
    for (ii = 0; ii < path_num; ii++)
    {
        req->path_list[ii] = path_list[ii];
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
                path->req_index = req->index;
            }
        }
    }
}

void init_rwa_db(RWA_DB *db, int linkNum, int pathNum, int reqGroupNum, int W)
{
    db->linkNum = linkNum;
    db->pathNum = pathNum;
    db->reqGroupNum = reqGroupNum;
    db->W = W;
    db->reqList = (RWA_REQ*)malloc(reqGroupNum * sizeof(RWA_REQ));
    db->pathList = (RWA_PATH*)malloc(pathNum * sizeof(RWA_PATH));
    db->linkList = (RWA_LINK*)malloc(linkNum * sizeof(RWA_LINK));
}

void del_db(RWA_DB *db)
{
    if (db != NULL)
    {
        del_link_list(db->linkList, db->linkNum);
        del_path_list(db->pathList, db->pathNum);
        del_req_list(db->reqList, db->reqGroupNum);
        free(db);
    }
}