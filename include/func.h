/**
 * @file func.h
 * @brief Function declaration
 * @author Liweining (liwn@cau.edu.cn)
 * @version 1.0
 * @date 2023-04-18
 *
 * @copyright Copyright (c) {2023}  Liu-Team
 *
 * @par change log:
 * <table>
 * <tr><th>Date       <th>Version <th>Author  <th>Description
 * <tr><td>2022-09-04 <td>1.0.0     <td>liwn     <td>get mrker information from
 * QMSim genotype file
 * </table>
 */

#ifdef linux
#include <sys/types.h>
#include <sys/wait.h>
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define FILE_BLOCK_SIZE 100 * 1024 * 1024

/* parse_opt 的参数 */
struct arguments
{
    int require_para;
    int onceWri;
    int nInd;
    int nMrk;
    char *mrkf;
    char *indexf;
    char *indIndexf;
    char *format;
    char *fid;
    char *out;
};

/* 根据索引输出基因型 */
int writeGT(int snp_total, int nsnp_sel, int onceWri, int nInd, int *mrk_index, char *format, char *fid, char *gtfile,
            char *outf, int *ind_index);
int countLines(char *file);
char *read_whole(char *file);
