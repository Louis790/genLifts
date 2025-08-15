#ifdef __cplusplus
extern "C" {
#endif

struct multigraph_result run_multigraph_gen(int argc, char** argv);
void reset();

struct multigraph_result {
    unsigned char* data;
    int graphs;
    int size;
};

#ifdef __cplusplus
}
#endif