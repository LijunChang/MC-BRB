#include "Utility.h"
#include "Graph.h"
#include "Timer.h"

void print_usage() {
	printf("Usage: [1]exe [2]alg [3]graph-dir [4 output]\n");
	printf("\talg: MC-DD, MC-EGO, MC-BRB, verify\n");
}

int main(int argc, char *argv[]) {
	if(argc < 3) {
		print_usage();
		return 0;
	}

#ifndef NDEBUG
	printf("**** MaxClique (Debug) build at %s %s ***\n", __TIME__, __DATE__);
#else
	printf("**** MaxClique (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

	printf("**** Alg: %s, Graph: %s ***\n", argv[1], argv[2]);

	printf("**** ");
#ifdef _KERNEL_
	printf("W/ kernelization, ");
#else
	printf("W/O kernelization, ");
#endif
#ifdef _RECOLOR_
	printf("W/ recolor, ");
#else
	printf("W/O recolor, ");
#endif
#ifdef _BRACH_ON_COLOR_
	printf("branch on color");
#endif
	printf(" ***\n");

	Graph *graph = new Graph(argv[2]);
	graph->read_graph_binary();
	//graph->read_graph(argv[2]);

#ifndef NDEBUG
	printf("\t*** Finished reading graph\n");
#endif

	if(strcmp(argv[1], "MC-DD") == 0) graph->heuristic_maximal_clique(1);
	else if(strcmp(argv[1], "MC-EGO") == 0) graph->maximum_clique_color(0,1);
	else if(strcmp(argv[1], "MC-BRB") == 0) graph->maximum_clique_color(1,1);
	else if(strcmp(argv[1], "verify") == 0) graph->verify_clique();
	else print_usage();

	if(argc >= 4&&strcmp(argv[3], "output") == 0) graph->output_clique();

	printf("\n");
	return 0;
}
