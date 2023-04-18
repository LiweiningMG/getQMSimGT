/**
 * @file func.c
 * @brief Function definition
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

#include "func.h"

/* 根据索引输出基因型 */
int writeGT(int snp_total, int nsnp_sel, int onceWri, int nInd, int *mrk_index, char *format, char *fid, char *gtfile,
            char *outf, int *ind_index)
{
    /* 基因型文件流 */
    FILE *QMSimGT = fopen(gtfile, "r");

    /* 输出文件 */
    FILE *output = fopen(outf, "w");

    /* arrays for storing genotypes */
    char *GT_all = (char *)calloc(snp_total * 4 + 100, sizeof(char));
    char *GT_mrk = (char *)calloc((nsnp_sel * 4 + 100) * onceWri, sizeof(char));

    /* 检查内存分配是否成功 */
    if (!GT_all || !GT_mrk)
    {
        printf("calloc memory fails!\n");
        return -1;
    }

    /* Gets the number of individuals that this thread needs to process */
    // int nInd_left = nInd;
    int snp_out = 0;
    int Indj_star_mrk = 0;
    int space = 0;
    int s = 1;
    int id_len = 0;
    char id[20] = "";
    char buf[100] = "";
    char header[100] = "";
    int headerL = 0;
    int nind = 0;

    /* 读取(跳过)标题行 */
    fgets(buf, 100, QMSimGT);

    /* 写出标题行 */
    if (strcmp(format, "ped") != 0)
    {
        fprintf(output, "%s", buf);
    }

    for (size_t indi = 0; indi < nInd; indi++)
    {
        /* 读取个体indi的基因型信息 */
        memset(GT_all, 0, sizeof(char) * (snp_total * 4 + 100));
        fgets(GT_all, snp_total * 4 + 100, QMSimGT);

        /* 判断该个体是否被选中输出基因型 */
        if (ind_index[indi] == 0)
        {
            continue;
        }

        /* 一组打包输出的个体计数 */
        nind++;

        /* 变量初始化 */
        space = 0;
        s = 1;
        snp_out = 0;

        /* 个体id */
        strcpy(id, strtok(GT_all, " "));
        id_len = strlen(id);

        /* 输出文件前n列 */
        memset(header, 0, sizeof(char) * 100);
        if (strcmp(format, "ped") == 0)
        {
            strcat(strcat(strcat(strcat(header, fid), " "), id), " 0 0 0 -9");
            headerL = id_len + strlen(fid) + 10;
        }
        else
        {
            strcpy(header, id);
            headerL = id_len;
        }

        /* 在输出文件中写入文件前n列 */
        for (int j = 0; j < headerL; j++)
        {
            GT_mrk[Indj_star_mrk + j] = header[j];
        }

        /* ID和基因型之间的空格数 */
        while (GT_all[id_len + s] == 32)
        {
            space++;
            s++;
        }

        /* 标记个体indi所在行的基因型信息起始位置 */
        Indj_star_mrk = headerL + Indj_star_mrk;

        /* 输出标记 */
        for (int j = 0; j < snp_total; j++)
        {
            /* 跳过不选择的标记 */
            if (mrk_index[j])
            {
                /* mrk genetype */
                GT_mrk[Indj_star_mrk + 4 * snp_out + 0] = 32; /* 32 is ASCII code for spaces */
                GT_mrk[Indj_star_mrk + 4 * snp_out + 1] = GT_all[4 * j + id_len + 1 + space];
                GT_mrk[Indj_star_mrk + 4 * snp_out + 2] = 32;
                GT_mrk[Indj_star_mrk + 4 * snp_out + 3] = GT_all[4 * j + id_len + 3 + space];

                snp_out++;
            }
        }

        /* 插入换行符 */
        GT_mrk[Indj_star_mrk + 4 * nsnp_sel] = 10;

        /* 下一个体信息的起始位置 */
        Indj_star_mrk = Indj_star_mrk + 4 * nsnp_sel + 1;

        /* 写出该组个体的基因型 */
        if (nind >= onceWri)
        {
            fprintf(output, "%s", GT_mrk);
            memset(GT_mrk, 0, sizeof(char) * (nsnp_sel * 4 + 100) * onceWri);

            /* 变量初始化 */
            nind = 0;
            Indj_star_mrk = 0;
        }
    }

    /* Free requested memory */
    free(GT_all);
    free(GT_mrk);
    fclose(output);
    fclose(QMSimGT);

    return 0;
}

char *read_whole(char *file)
{
    /* Read in binary */
    FILE *pFile = fopen(file, "rb");
    if (pFile == NULL)
    {
        printf("Can not load file: %s\n", file);
        exit(1);
    }

    /* get file size */
    fseek(pFile, 0, SEEK_END);
    long lSize = ftell(pFile);
    rewind(pFile);

    /* allocate memory */
    char *buffer = (char *)calloc(lSize, sizeof(char));
    if (buffer == NULL)
    {
        fputs("Memory allocate error.", stderr);
        exit(2);
    }

    /* read whole file to buffer */
    size_t result = fread(buffer, 1, lSize, pFile);
    if (result != lSize)
    {
        fputs("fread error", stderr);
        exit(3);
    }

    fclose(pFile);
    return buffer;
}

#ifdef _WIN32
int countLines(char *file)
{
    FILE *fp;
    int rows = 0;
    char *BufContent = (char *)calloc(FILE_BLOCK_SIZE, sizeof(char));
    size_t BufContentSz;

    if ((fp = fopen(file, "rb")) == NULL)
    {
        printf("Can not load file: %s\n", file);
        exit(EXIT_FAILURE);
    }
    if (BufContent == NULL)
    {
        printf("The memory of size %d cannot be allocated in the count function.\n", FILE_BLOCK_SIZE);
        fclose(fp);
        return -2;
    }
    else
    {
        while ((BufContentSz = fread(BufContent, sizeof(unsigned char), FILE_BLOCK_SIZE, fp)) > 0)
        {
            for (int i = 0; i < BufContentSz; i++)
            {
                if (BufContent[i] == '\n')
                {
                    rows++;
                }
            }
        }
        free(BufContent);
    }

    free(BufContent);
    return rows;
}
#endif

#ifdef linux
int countLines(char *file)
{
    int rows = 0;
    char cmdstring[210];
    memset(cmdstring, 0, sizeof(cmdstring));
    char buff[1024];
    memset(buff, 0, sizeof(buff));
    FILE *fstream = NULL;

    /*准备输出文件行数的bash命令*/
    strcpy(buff, file);
    strcat(buff, " | awk '{print $1}'");
    strcat(cmdstring, "wc -l ");
    strcat(cmdstring, buff);

    if ((fstream = popen(cmdstring, "r")) == NULL)
    {
        fprintf(stderr, "execute command failed: %s", strerror(errno));
        return -1;
    }

    if (fgets(buff, sizeof(buff), fstream) != NULL)
    {
        rows = atoi(buff);
    }
    else
    {
        rows = -1;
    }

    pclose(fstream);
    return rows;
}
#endif
