#include <stdio.h>
#include <stdlib.h>
#include <argp.h>

#include "func.h"

/* 程序说明 */
const char *argp_program_version = "getQMSimPed 1.0";
const char *argp_program_bug_address = "<liwn@cau.edu.cn>";
static char doc[] = "QMSim_selected -- Output the markers selected";

/* parse function */
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
  /* 从 argp_parse 获取参数，该参数为 arguments 结构体 */
  struct arguments *arguments = state->input;

switch (key) {
  case 'm':
    arguments->mrkf = arg;
    arguments->require_para++;
    break;
  case 'i':
    arguments->indexf = arg;
    arguments->require_para++;
    break;
  case 'f':
    arguments->format = arg;
    break;
  case 'd':
    arguments->fid = arg;
    break;
  case 'w':
    arguments->onceWri = atoi(arg);
    break;
  case 'n':
    arguments->nInd = atoi(arg);
    break;
  case 'N':
    arguments->nMrk = atoi(arg);
    break;
  case 'I':
    arguments->indIndexf = arg;
    break;
  case 'o':
    arguments->out = arg;
    break;
  case ARGP_KEY_END: {
    if (arguments->require_para < 2) {
      printf("Missing required parameters [%d].\n", arguments->require_para);
      argp_state_help(state, state->err_stream, ARGP_HELP_STD_HELP);
      exit(-1);
    }
  } break;
  default:
    return ARGP_ERR_UNKNOWN;
  }

  return 0;
}

/* Parameter description */
static struct argp_option options[] = {
    {0, 0, 0, 0, "Required parameters:", 1},
    {"mrkf", 'm', "FILE", 0, "**_mrk_001.txt"},
    {"indexf", 'i', "FILE", 0, "**.txt"},
    {0, 0, 0, 0, "Optional parameters:", 2},
    {"indIndexf", 'I', "FILE", 0, "Location of the selected individual [NULL]"},
    {"nInd", 'n', "INT", 0, "Number of individuals [0]"},
    {"nMrk", 'N', "INT", 0, "Number of markers [0]"},
    {"onceWri", 'w', "INT", 0, "Number of individuals written out to the file once [50]"},
    {"format", 'f', "CHR", 0, "mrk/ped [ped]"},
    {"fid", 'd', "CHR", 0, "family id in ped file [breedA]"},
    {"out", 'o', "FILE", 0, "[QMSim_selected.ped]"},
    {0}};

/* argp parser, parameters of argp_parse */
static struct argp argp = {options, parse_opt, 0, doc};

int main(int argc, char *argv[])
{
    /* ----------- 命令行参数处理 ----------- */
    struct arguments para;

    /* 可选参数默认值 */
    para.nInd = 0;
    para.nMrk = 0;
    para.onceWri = 50;
    para.indIndexf = NULL;
    para.fid = (char *)"breedA";
    para.format = (char *)"ped";
    para.out = (char *)"QMSim_selected.ped";

    /* 解析命令行参数  */
    para.require_para = 0; /* 必要参数个数 */
    argp_parse(&argp, argc, argv, ARGP_IN_ORDER, 0, &para);

    /* 基因型文件中个体数 */
    if (para.nInd == 0)
    {
        printf("Counting the number of individuals in mrkf...\n");
        para.nInd = countLines(para.mrkf) - 1;
        printf("number of individuals in mrkf: %d\n", para.nInd);
    }

    /* 标记数目 */
    if (para.nMrk == 0)
    {
        printf("Counting the number of markers...\n");
        para.nMrk = countLines(para.indexf);
        printf("number of markers in mrkf: %d\n", para.nMrk);
    }

    /* 读取标记选择索引文件(0/1) */
    char *index_chr = read_whole(para.indexf);
    int *mrk_index = (int *)calloc(para.nMrk, sizeof(int));
    int nmrk_sel = 0;
    for (int i = 0; i < para.nMrk; i++)
    {
        /* 选择的SNP标记数 */
        if (strcmp(index_chr + 2 * i, "1") == 10)
        {
            nmrk_sel++;
            mrk_index[i] = 1;
        }
    }
    printf("number of marker selected: %d\n", nmrk_sel);

    /* 读取个体选择索引文件(0/1) */
    int *ind_index = (int *)calloc(para.nInd, sizeof(int));
    int nind_sel = 0;
    if (strlen(para.indIndexf) == 0)
    {
        /* 全选 */
        nind_sel = para.nInd;
        for (int i = 0; i < para.nInd; i++)
        {
            ind_index[i] = 1;
        }
    }
    else
    {
        // 检查选择个体的索引数是不是和个体数一致
        if (countLines(para.indIndexf) != para.nInd)
        {
            printf("\nError: The number of individuals in indIndexf is not equal to that in mrkf!\n");
            exit(1);
        }

        char *ind_index_chr = read_whole(para.indIndexf);

        for (int i = 0; i < para.nInd; i++)
        {
            /* 选择的SNP标记数 */
            if (strcmp(ind_index_chr + 2 * i, "1") == 10)
            {
                nind_sel++;
                ind_index[i] = 1;
            }
        }
    }

    printf("number of individuals selected: %d\n", nind_sel);

    /* 输出基因型文件 */
    printf("writing out the file with %d individuals at a time\n", para.onceWri);
    writeGT(para.nMrk, nmrk_sel, para.onceWri, para.nInd, mrk_index, para.format, para.fid, para.mrkf, para.out, ind_index);

    free(mrk_index);
    free(ind_index);

    return 0;
}
