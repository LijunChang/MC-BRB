#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "LinearHeap.h"

class Graph {
private:
	std::string dir; //input graph directory
	ui n; //number of nodes of the graph
	ept m; //number of edges of the graph

	ept *pstart; //offset of neighbors of nodes
	ept *pend; //used in search
	ui *edges; //adjacent ids of edges

	ui *tmp_vs;
	ui *tmp_color;

	ept *pstart_o; //oriented graph
	ui *edges_o; //oriented graph

	unsigned char *matrix; //adjacency matrix of an induced subgraph

	ui *head; // size of max_core, used for sorting vertices w.r.t. #colors
	ui *next; // size of max_core, used for sorting vertices w.r.t. #colors
	ui *id;

	ui max_core;
	ui max_depth;
	long long branches;

	ui *degree;
	ui *rid;
	char *vis;

	ui *current_clique;
	std::vector<ui> vs_buf;
	std::vector<ui> color_buf;

	ui *mapping;
	ui mapping_n;
	ui matrix_len;

	std::vector<ui> max_clique;
	std::vector<ui> contractions;
	std::vector<std::pair<ui, ui> > changes;
	std::vector<ui> del, degree_one, degree_two, degree_three;

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph_binary() ;
	void read_graph(const char *input_file) ;
	void maximum_clique_color(char exact, char initial_by_MC_EDGE) ;
	void heuristic_maximal_clique(char opt) ;
	void output_clique() ;
	void verify_clique() ;

private:
	ui coloring_adj_list(const ui *vs, const ui vs_size, const ui original_size, ui *color, char *vis, const ui start_idx, const ui start_color) ;
	ui degeneracy_ordering_and_coloring_adj(ui *vs, ui *color, const ui vs_size, char *vis, ui *degree, ui *label, const ui K, const ui threshold) ;
	ui degeneracy_maximal_clique_adjacency_list(ui *peel_sequence, ui *core, ui *color, char opt, char greedy_extend = 0) ;
	ui degeneracy_maximal_clique_adj_list(const ui current_clique_size, ui *vs, const ui vs_size, char *vis, ui *degree, ListLinearHeap *heap) ;
	ui ego_degen(const ui *peel_sequence, const ui *core, const ui *color, ui *local_UBs, const ui UB) ;

	// build the oriented graph and store in pstart_o and edges_o, mapping of ids is stored in out_mapping
	// also reduces the sizes of peel_sequence, core, and color
	void shrink_graph(ui *&peel_sequence, ui *&core, ui *&color, ui *&out_mapping, ept *&pstart_o, ui *&edges_o) ;

	void search_oriented(const ui *peel_sequence, const ui *core, const ui *color, const ui *local_UBs) ;

	void recursive_search_clique_color_with_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size, char kernel) ;
	void recursive_search_clique_color_without_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size) ;

	void heuristic_max_clique_max_degree(ui processed_threshold) ;

	ui color_bound(const ui *vs, const ui vs_size, const ui *color, char *vis) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) assert(!vis[color[vs[i]]]);
#endif
		ui color_bound = 0;
		for(ui i = 0;i < vs_size;i ++) if(!vis[color[vs[i]]]) {
			vis[color[vs[i]]] = 1;
			++ color_bound;
		}
		for(ui i = 0;i < vs_size;i ++) vis[color[vs[i]]] = 0;
		return color_bound;
	}

	ui color_bound(const ui *vs, const ui vs_size, const ui *mapping, const ui *color, char *vis) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) assert(!vis[color[mapping[vs[i]]]]);
#endif
		ui color_bound = 0;
		for(ui i = 0;i < vs_size;i ++) if(!vis[color[mapping[vs[i]]]]) {
			vis[color[mapping[vs[i]]]] = 1;
			++ color_bound;
		}
		for(ui i = 0;i < vs_size;i ++) vis[color[mapping[vs[i]]]] = 0;
		return color_bound;
	}

	//coloring a graph that is represented by matrix
	//return the number of colors used
	ui coloring_matrix(const ui *vs, const ui vs_size, ui *color, char *vis, const ui start_idx, const ui start_color) {
		assert(start_color > 0&&start_color <= vs_size);
		for(ui i = vs_size - start_color;i < vs_size;i ++) color[vs[i]] = vs_size - i - 1;
		ui max_color = start_color-1;

		for(ui i = vs_size - start_color;i > start_idx;i --) {
			ui u = vs[i-1];
			unsigned char *t_matrix = matrix + u*matrix_len;
			for(ui j = i;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) vis[color[vs[j]]] = 1;
			for(ui j = 0;;j ++) if(!vis[j]) {
				color[u] = j;
				if(j > max_color) max_color = j;
				break;
			}
			for(ui j = i;j < vs_size;j ++) vis[color[vs[j]]] = 0;
		}

		return max_color + 1;
	}

	//coloring a graph that is represented by matrix, aiming to minimize the number of vertices with color >= threshold
	//return the number of colors used
	ui coloring_matrix_advanced(const ui *vs, const ui vs_size, ui *color, const ui start_color, const ui threshold) {
		assert(rid != nullptr&&head != nullptr); // rid is used to temporarily store the color
		for(ui i = 0;i < vs_size;i ++) head[i] = n;

		assert(start_color > 0&&start_color <= vs_size);
		for(ui i = vs_size - start_color;i < vs_size;i ++) {
			ui c = vs_size - i - 1;
			rid[vs[i]] = c;
			id[i] = vs[i];
			next[i] = head[c];
			head[c] = i;
		}

		for(ui i = vs_size - start_color;i > 0;) {
			ui u = vs[-- i];
			unsigned char *t_matrix = matrix + u*matrix_len;
			rid[u] = n;
			for(ui j = 0;j < threshold;j ++) {
				char ok = 1;
				for(ui k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
					ok = 0;
					break;
				}
				if(ok) {
					rid[u] = j;
					id[i] = vs[i];
					next[i] = head[j];
					head[j] = i;
					break;
				}
			}
#ifndef _RECOLOR_
			continue;
#endif
			if(rid[u] < threshold) continue;

			for(ui j = 0;j < threshold;j ++) {
				ui cnt = 0, idx;
				for(ui k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
					++ cnt;
					if(cnt == 1) idx = k;
					else break;
				}
				assert(cnt > 0);
				if(cnt != 1) continue;
				unsigned char *tt_matrix = matrix + id[idx]*matrix_len;
				for(ui ii = threshold;ii > 0;) {
					-- ii;
					if(ii == j) continue;

					char ok = 1;
					for(ui k = head[ii];k != n;k = next[k]) if(test_bit(tt_matrix, id[k])) {
						ok = 0;
						break;
					}
					if(ok) {
						rid[id[idx]] = ii;
						id[i] = id[idx];
						next[i] = head[ii];
						head[ii] = i;

						rid[u] = j;
						id[idx] = vs[i];
						break;
					}
				}
				if(rid[u] < threshold) break;
			}
			if(rid[u] < threshold) continue;

			for(ui j = threshold;;j ++) {
				char ok = 1;
				for(ui k = head[j];k != n;k = next[k]) if(test_bit(t_matrix, id[k])) {
					ok = 0;
					break;
				}
				if(ok) {
					rid[u] = j;
					id[i] = vs[i];
					next[i] = head[j];
					head[j] = i;
					break;
				}
			}
		}

		ui max_color = 0;
		for(ui i = 0;i < vs_size;i ++) {
			color[i] = rid[vs[i]];
			if(color[i] > max_color) max_color = color[i];
		}

#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) {
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) assert(color[i] != color[j]);
		}
#endif
		return max_color+1;
	}

	// construct induced subgraph from pstart_o and edges_o
	void construct_induced_subgraph(const ui *vs, const ui vs_size, char *vis, ui *degree, const ept *pstart_o, const ui *edges_o) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) assert(vis[vs[i]]);
#endif
		for(ui j = 0;j < vs_size;j ++) degree[vs[j]] = 0;
		for(ui j = 0;j < vs_size;j ++) {
			ui v = vs[j];
			for(ui k = pstart_o[v];k < pstart_o[v+1];k ++) {
				if(vis[edges_o[k]]) {
					++ degree[v];
					++ degree[edges_o[k]];
				}
			}
		}

		pstart[vs[0]] = 0;
		for(ui j = 1;j < vs_size;j ++) pstart[vs[j]] = pstart[vs[j-1]] + degree[vs[j-1]];
		for(ui j = 0;j < vs_size;j ++) pend[vs[j]] = pstart[vs[j]];
		for(ui j = 0;j < vs_size;j ++) {
			ui v = vs[j];
			for(ui k = pstart_o[v];k < pstart_o[v+1];k ++) {
				ui w = edges_o[k];
				if(vis[w]) {
					edges[pend[v] ++] = w;
					edges[pend[w] ++] = v;
				}
			}
		}
	}

	// construct matrix from pstart_o and edges_o for vertices in vs
	void construct_matrix(ui *vs, const ui vs_size, ui *mapping, char *vis, ui *rdegree, const ept *pstart_o, const ui *edges_o) {
		assert(rid != nullptr&&mapping != nullptr&&matrix != nullptr);
		assert(vs_size <= max_core);

		mapping_n = vs_size;
#ifdef _BITSET_
		matrix_len = (vs_size+7)/8;
#else
		matrix_len = vs_size;
#endif
		for(ui j = 0;j < vs_size;j ++) {
			vis[vs[j]] = 1;
			rid[vs[j]] = j;
			mapping[j] = vs[j];
			vs[j] = j;
		}

		memset(matrix, 0, sizeof(unsigned char)*mapping_n*matrix_len);
		for(ui j = 0;j < mapping_n;j ++) rdegree[j] = mapping_n - 1;
		for(ui j = 0;j < mapping_n;j ++) for(ui k = pstart_o[mapping[j]];k < pstart_o[mapping[j]+1];k ++) if(vis[edges_o[k]]) {
			ui v = rid[edges_o[k]];
			assert(v >= 0&&v < mapping_n);
			set_bit(matrix + j*matrix_len, v);
			set_bit(matrix + v*matrix_len, j);
			-- rdegree[j]; -- rdegree[v];
		}
		for(ui j = 0;j < mapping_n;j ++) set_bit(matrix + j*matrix_len, j);

		for(ui j = 0;j < mapping_n;j ++) vis[mapping[j]] = 0;
	}

	ui degeneracy_maximal_clique_matrix(ui current_clique_size, ui *vs, const ui vs_size, ui *degree, char heuristic_gen, char print) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) {
			ui d = 0;
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = 0;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) ++ d;
			assert(test_bit(t_matrix, vs[i]));
			assert(degree[vs[i]]+1 == d);
		}
#endif

		ui start_color = 0;
		for(ui j = 0;j < vs_size;j ++) {
			ui min_idx = j;
			for(ui k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
			if(min_idx != j) std::swap(vs[min_idx], vs[j]);
			ui v = vs[j];
			if(degree[v] + 1 + j == vs_size) {
				for(ui k = j+1;k < vs_size;k ++) degree[vs[k]] = vs_size-1-k;
				start_color = degree[v];
				assert(start_color > 0&&start_color < vs_size);
				if(heuristic_gen&&degree[v] + 1 + current_clique_size > max_clique.size()) {
					for(ui k = j;k < vs_size;k ++) current_clique[current_clique_size ++] = vs[k];
					store_a_larger_clique(current_clique_size, "degen", print);
				}
				break;
			}
			unsigned char *t_matrix = matrix + v*matrix_len;
			for(ui k = j+1;k < vs_size;k ++) if(test_bit(t_matrix, vs[k])) -- degree[vs[k]];
		}
		return start_color;
	}

	void degree_one_two_reduction_with_folding_matrix(ui &current_clique_size, ui *vs, ui &vs_size, ui *rdegree, std::vector<ui> &contractions) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) {
			ui rd = 0;
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) ++ rd;
			assert(test_bit(t_matrix, vs[i]));
			assert(rdegree[vs[i]] == rd);
		}
#endif

		char changed = 1;
		while(changed) {
			changed = 0;
			for(ui i = 0;i < vs_size;) {
				ui rd = rdegree[vs[i]];
				if(rd == 0) {
					changed = 1;
					current_clique[current_clique_size ++] = vs[i];
					vs[i] = vs[-- vs_size];
				}
				else if(vs_size - rd + current_clique_size <= max_clique.size()) {
					changed = 1;
					remove_vertex_idx(vs, i, vs_size, rdegree);
				}
				else if(rd == 1) {
					changed = 1;
					current_clique[current_clique_size ++] = vs[i];
					unsigned char *t_matrix = matrix + vs[i]*matrix_len;
					vs[i] = vs[-- vs_size];
					for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
						remove_vertex_idx(vs, j, vs_size, rdegree);
						break;
					}
				}
				else if(rd == 2) {
					changed = 1;
					ui u = vs[i];
					unsigned char *t_matrix = matrix + u*matrix_len;
					vs[i] = vs[-- vs_size];
					ui idx1 = n, idx2 = n;
					for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
						if(idx1 == n) idx1 = j;
						else {
							idx2 = j;
							break;
						}
					}
					assert(idx1 != n&&idx2 != n&&idx1 < idx2);

					if(!test_bit(matrix+vs[idx1]*matrix_len, vs[idx2])) {
						current_clique[current_clique_size ++] = u;
						remove_vertex_idx(vs, idx2, vs_size, rdegree);
						remove_vertex_idx(vs, idx1, vs_size, rdegree);
					}
					else {
						current_clique[current_clique_size ++] = n;
						-- rdegree[vs[idx1]];
						contractions.pb(vs[idx1]); contractions.pb(vs[idx2]); contractions.pb(u); contractions.pb(1);
						unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
						unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
						vs[idx2] = vs[-- vs_size];

						for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix2, vs[j])) {
							if(test_bit(t_matrix1, vs[j])) {
								reverse_bit(t_matrix1, vs[j]);
								reverse_bit(matrix + vs[j]*matrix_len, vs[idx1]);
								++ rdegree[vs[idx1]];
							}
							else -- rdegree[vs[j]];
						}
					}
				}
				else ++ i;
			}
		}

#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) {
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			ui rd = 0;
			for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) ++ rd;
			assert(rd == rdegree[vs[i]]);
			assert(vs_size - 1 - rd + current_clique_size >= max_clique.size());
			assert(rd >= 2);
		}
#endif

		if(current_clique_size > max_clique.size()) store_a_larger_clique(current_clique_size, "degree_one_two_with_folding", 0);
	}

	void put_into_one_vector(ui &current_clique_size, ui &i, ui *vs, ui &vs_size, bool check, const ui rd, ui *rid) {
		if(vs_size - rd + current_clique_size <= max_clique.size()) {
			del.pb(vs[i]);
			//printf("insert %u\n", vs[i]);
		}
		else if(check&&rd <= 3) {
			switch(rd) {
			case 0: current_clique[current_clique_size ++] = vs[i];
					vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
					-- i;
					break;
			case 1: degree_one.pb(vs[i]);
					break;
			case 2: degree_two.pb(vs[i]);
					break;
			case 3: degree_three.pb(vs[i]);
			}
		}
	}

	void put_into_one_vector_eq(ui &current_clique_size, ui &i, ui *vs, ui &vs_size, bool check, const ui rd, ui *rid) {
		if(check) {
			if(rd <= 3&&vs_size - rd + current_clique_size > max_clique.size()) {
				switch(rd) {
				case 0: current_clique[current_clique_size ++] = vs[i];
						//printf("added %u to maximum clique\n", vs[i]);
						vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
						-- i;
						break;
				case 1: degree_one.pb(vs[i]);
					break;
				case 2: degree_two.pb(vs[i]);
					break;
				case 3: degree_three.pb(vs[i]);
				}
			}
		}
		else {
			if(vs_size - rd + current_clique_size == max_clique.size()) {
				del.pb(vs[i]);
				//printf("insert %u\n", vs[i]);
			}
		}
	}

	void degree_one_two_three_reduction_with_folding_matrix(ui &current_clique_size, ui *vs, ui &vs_size, ui *rdegree, ui *rid) {
		assert(del.empty()&&degree_one.empty()&&degree_two.empty()&&degree_three.empty());

		for(ui i = 0;i < vs_size;i ++) {
			rid[vs[i]] = i;
			put_into_one_vector(current_clique_size, i, vs, vs_size, true, rdegree[vs[i]], rid);
		}

		while(!del.empty()||!degree_one.empty()||!degree_two.empty()||!degree_three.empty()) {
			while(!del.empty()) {
				ui u = del.back(); del.pop_back();
				rdegree[u] = 0;
				//printf("u: %u, rid[u]: %u, vs[rid[u]]: %u, vs_size: %u\n", u, rid[u], vs[rid[u]], vs_size);
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				ui idx = rid[u];
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;
				unsigned char *t_matrix = matrix + u*matrix_len;
				for(ui i = 0;i < vs_size;i ++) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;
					if(!test_bit(t_matrix, vs[i])) -- rd;
					put_into_one_vector_eq(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}

#ifndef NDEBUG
				for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}
			while(del.empty()&&!degree_one.empty()) {
				ui u = degree_one.back(); degree_one.pop_back();
				if(rdegree[u] != 1) continue;

				current_clique[current_clique_size ++] = u;
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				ui idx = rid[u]; rdegree[u] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

				unsigned char *t_matrix = matrix + u*matrix_len; idx = n;
				for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					idx = j;
					break;
				}
				assert(idx != n);

				t_matrix = matrix + vs[idx]*matrix_len;
				rdegree[vs[idx]] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;
				for(ui i = 0;i < vs_size;i ++) {
					ui &rd = rdegree[vs[i]];
					ui old_rd = rd;
					if(!test_bit(t_matrix, vs[i])) -- rd;
					put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
				}

#ifndef NDEBUG
				for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}
			while(del.empty()&&degree_one.empty()&&!degree_two.empty()) {
				ui u = degree_two.back(); degree_two.pop_back();
				if(rdegree[u] != 2) continue;

				unsigned char *t_matrix = matrix + u*matrix_len;
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				ui idx = rid[u]; rdegree[u] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

				ui idx1 = n, idx2 = n;
				for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					if(idx1 == n) idx1 = j;
					else {
						idx2 = j;
						break;
					}
				}
				assert(idx1 != n&&idx2 != n&&idx1 < idx2);

				unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
				unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
				if(!test_bit(t_matrix1, vs[idx2])) {
					current_clique[current_clique_size ++] = u;
					rdegree[vs[idx2]] = 0;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					rdegree[vs[idx1]] = 0;
					vs[idx1] = vs[-- vs_size]; rid[vs[idx1]] = idx1;

					for(ui i = 0;i < vs_size;i ++) {
						ui &rd = rdegree[vs[i]];
						ui old_rd = rd;
						if(!test_bit(t_matrix1, vs[i])) -- rd;
						if(!test_bit(t_matrix2, vs[i])) -- rd;

						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
				}
				else {
					current_clique[current_clique_size ++] = n;
					ui v = vs[idx1], old_degree_v = rdegree[v];
					-- rdegree[v];
					contractions.pb(v); contractions.pb(vs[idx2]); contractions.pb(u); contractions.pb(1); //type-1 contraction

					rdegree[vs[idx2]] = 0;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;

					for(ui i = 0;i < vs_size;i ++) if(vs[i] != v) {
						ui &rd = rdegree[vs[i]];
						ui old_rd = rd;

						if(!test_bit(t_matrix2, vs[i])) {
							if(test_bit(t_matrix1, vs[i])) {
								reverse_bit(t_matrix1, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, v);
								++ rdegree[v];
								changes.pb(std::make_pair(v, vs[i]));
							}
							else -- rd;
						}
						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
				}

#ifndef NDEBUG
				for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}

			while(del.empty()&&degree_one.empty()&&degree_two.empty()&&!degree_three.empty()) {
				ui u = degree_three.back(); degree_three.pop_back();
				if(rdegree[u] != 3) continue;

				unsigned char *t_matrix = matrix + u*matrix_len;
				assert(rid[u] < vs_size&&vs[rid[u]] == u);
				ui idx = rid[u];
				//std::swap(vs[idx], vs[vs_size-1]); rid[vs[idx]] = idx; rid[vs[vs_size-1]] = vs_size - 1;
				rdegree[u] = 0;
				vs[idx] = vs[-- vs_size]; rid[vs[idx]] = idx;

				ui idx1 = n, idx2 = n, idx3 = n;
				for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) {
					if(idx1 == n) idx1 = j;
					else if(idx2 == n) idx2 = j;
					else {
						idx3 = j;
						break;
					}
				}
				assert(idx1 != n&&idx2 != n&&idx3 != n&&idx1 < idx2&&idx2 < idx3);

				unsigned char *t_matrix1 = matrix + vs[idx1]*matrix_len;
				unsigned char *t_matrix2 = matrix + vs[idx2]*matrix_len;
				unsigned char *t_matrix3 = matrix + vs[idx3]*matrix_len;
				char connected_1 = 0, connected_2 = 0, connected_3 = 0;
				if(test_bit(t_matrix1, vs[idx2])) connected_1 = 1;
				if(test_bit(t_matrix1, vs[idx3])) connected_2 = 1;
				if(test_bit(t_matrix2, vs[idx3])) connected_3 = 1;
				char total_connected = connected_1 + connected_2 + connected_3;
				if(total_connected == 0) { //isolation reduction
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = u;
					rdegree[vs[idx3]] = 0;
					vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
					rdegree[vs[idx2]] = 0;
					vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					rdegree[vs[idx1]] = 0;
					vs[idx1] = vs[-- vs_size]; rid[vs[idx1]] = idx1;

					for(ui i = 0;i < vs_size;i ++) {
						ui &rd = rdegree[vs[i]];
						ui old_rd = rd;
						if(!test_bit(t_matrix1, vs[i])) -- rd;
						if(!test_bit(t_matrix2, vs[i])) -- rd;
						if(!test_bit(t_matrix3, vs[i])) -- rd;

						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
				}
				else if(total_connected == 1) {
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = n;
					//the following ensures that vs[idx2] is the vertex to be deleted
					if(connected_3) {
						std::swap(vs[idx1], vs[idx2]);
						std::swap(t_matrix1, t_matrix2);
						rid[vs[idx1]] = idx1;
					}
					else if(connected_1) {
						std::swap(vs[idx2], vs[idx3]);
						std::swap(t_matrix2, t_matrix3);
					}

					ui v = vs[idx1], old_degree_v = rdegree[v];
					rdegree[v] -= 2;
					contractions.pb(v); contractions.pb(vs[idx3]); contractions.pb(u); contractions.pb(1); //type-1 contraction

					rdegree[vs[idx3]] = 0; rdegree[vs[idx2]] = 0;
					assert(idx2 != idx3);
					if(idx2 < idx3) {
						vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
						vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
					}
					else {
						vs[idx2] = vs[-- vs_size]; rid[vs[idx2]] = idx2;
						vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;
					}

					for(ui i = 0;i < vs_size;i ++) if(vs[i] != v) {
						ui &rd = rdegree[vs[i]];
						ui old_rd = rd;

						if(!test_bit(t_matrix2, vs[i])) -- rd;

						if(!test_bit(t_matrix3, vs[i])) {
							if(test_bit(t_matrix1, vs[i])) {
								reverse_bit(t_matrix1, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, v);
								changes.pb(std::make_pair(v, vs[i]));
								++ rdegree[v];
							}
							else -- rd;
						}
						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
				}
				else if(total_connected == 2) {
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = n;
					//the following ensures that vs[idx2] is the vertex to be deleted
					if(!connected_3) {
						std::swap(vs[idx1], vs[idx3]);
						std::swap(t_matrix1, t_matrix3);
						rid[vs[idx1]] = idx1;
					}
					else if(!connected_2) {
						std::swap(vs[idx2], vs[idx3]);
						std::swap(t_matrix2, t_matrix3);
						rid[vs[idx2]] = idx2;
					}

					ui v = vs[idx1], old_degree_v = rdegree[v];
					ui w = vs[idx2], old_degree_w = rdegree[w];
					-- rdegree[v]; -- rdegree[w];
					contractions.pb(v); contractions.pb(w); contractions.pb(vs[idx3]); contractions.pb(u); contractions.pb(2); //type-2 contraction

					rdegree[vs[idx3]] = 0;
					vs[idx3] = vs[-- vs_size]; rid[vs[idx3]] = idx3;

					for(ui i = 0;i < vs_size;i ++) if(vs[i] != v&&vs[i] != w) {
						ui &rd = rdegree[vs[i]];
						ui old_rd = rd;

						if(!test_bit(t_matrix3, vs[i])) {
							-- rd;
							if(test_bit(t_matrix1, vs[i])) {
								reverse_bit(t_matrix1, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, v);
								changes.pb(std::make_pair(v, vs[i]));
								++ rdegree[v];
								++ rd;
							}

							if(test_bit(t_matrix2, vs[i])) {
								reverse_bit(t_matrix2, vs[i]);
								reverse_bit(matrix + vs[i]*matrix_len, w);
								changes.pb(std::make_pair(w, vs[i]));
								++ rdegree[w];
								++ rd;
							}
						}
						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v], vs, vs_size, rdegree[v] != old_degree_v, rdegree[v], rid);
					put_into_one_vector(current_clique_size, rid[w], vs, vs_size, rdegree[w] != old_degree_w, rdegree[w], rid);
				}
				else {
					assert(total_connected == 3);
					//rdegree[u] = 0; -- vs_size;

					current_clique[current_clique_size ++] = n;

					ui v1 = vs[idx1], old_degree_v1 = rdegree[v1];
					ui v2 = vs[idx2], old_degree_v2 = rdegree[v2];
					ui v3 = vs[idx3], old_degree_v3 = rdegree[v3];
					-- rdegree[v3];
					reverse_bit(t_matrix1, v2); reverse_bit(t_matrix2, v1);
					changes.pb(std::make_pair(v1, v2));
					contractions.pb(v1); contractions.pb(v2); contractions.pb(v3); contractions.pb(u); contractions.pb(3); //type-3 contraction

					for(ui i = 0;i < vs_size;i ++) if(vs[i] != v1&&vs[i] != v2&&vs[i] != v3) {
						ui &rd = rdegree[vs[i]];
						ui old_rd = rd;

						char conn1 = test_bit(t_matrix1, vs[i]);
						char conn2 = test_bit(t_matrix2, vs[i]);
						char conn3 = test_bit(t_matrix3, vs[i]);

						if(conn1&&!conn2) {
							reverse_bit(t_matrix1, vs[i]);
							reverse_bit(matrix+vs[i]*matrix_len, v1);
							changes.pb(std::make_pair(v1, vs[i]));
							++ rd;
							++ rdegree[v1];
						}

						if(conn2&&!conn3) {
							reverse_bit(t_matrix2, vs[i]);
							reverse_bit(matrix+vs[i]*matrix_len, v2);
							changes.pb(std::make_pair(v2, vs[i]));
							++ rd;
							++ rdegree[v2];
						}

						if(conn3&&!conn1) {
							reverse_bit(t_matrix3, vs[i]);
							reverse_bit(matrix+vs[i]*matrix_len, v3);
							changes.pb(std::make_pair(v3, vs[i]));
							++ rd;
							++ rdegree[v3];
						}

						put_into_one_vector(current_clique_size, i, vs, vs_size, old_rd != rd, rd, rid);
					}
					put_into_one_vector(current_clique_size, rid[v1], vs, vs_size, rdegree[v1] != old_degree_v1, rdegree[v1], rid);
					put_into_one_vector(current_clique_size, rid[v2], vs, vs_size, rdegree[v2] != old_degree_v2, rdegree[v2], rid);
					put_into_one_vector(current_clique_size, rid[v3], vs, vs_size, rdegree[v3] != old_degree_v3, rdegree[v3], rid);
				}

#ifndef NDEBUG
				for(ui i = 0;i < vs_size;i ++) assert(rid[vs[i]] == i);
#endif
			}
		}

#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) assert(rdegree[vs[i]] > 3&&vs_size - rdegree[vs[i]] + current_clique_size > max_clique.size());
#endif
	}

	void get_higher_neighbors(const ui u, ui &vs_size, std::vector<ui> &vs_buf, std::vector<ui> &color_buf, const ept *pstart_o, const ui *edges_o) {
		vs_size = 0;
		for(ui j = pstart_o[u];j < pstart_o[u+1];j ++) {
			if(vs_buf.size() == vs_size) {
				vs_buf.pb(edges_o[j]);
				color_buf.pb(0);
			}
			else vs_buf[vs_size] = edges_o[j];
			++ vs_size;
		}
	}

	//greedily enlarge max_clique by including u if feasible
	char greedy_extend(const ui u, const char *vis, char print) {
		for(auto v: max_clique) if(!vis[v]) return 0;
		max_clique.pb(u);
		if(print) printf("greedy_extend finds clique of size: %lu\n", max_clique.size());
		return 1;
	}

	void kcore_reduction(ui *vs, ui &vs_size, char *vis, ui *degree, const ui K, ui *queue) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) assert(vis[vs[i]]);
#endif

		ui queue_n = 0;
		for(ui i = 0;i < vs_size;i ++) if(degree[vs[i]] < K) {
			vis[vs[i]] = 0;
			queue[queue_n ++] = vs[i];
		}
		for(ui i = 0;i < queue_n;i ++) {
			ui u = queue[i];
			for(ui j = pstart[u];j < pend[u];j ++) if(vis[edges[j]]) {
				if((-- degree[edges[j]]) == K-1) {
					vis[edges[j]] = 0;
					queue[queue_n ++] = edges[j];
				}
			}
		}
		ui new_vs_size = 0;
		for(ui i = 0;i < vs_size;i ++) if(vis[vs[i]]) {
			if(new_vs_size != i) std::swap(vs[new_vs_size], vs[i]);
			++ new_vs_size;
		}
		vs_size = new_vs_size;
	}


	char kernelization_color(ui &clique_size, ui *vs, ui &vs_size, ui *color) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size-1;i ++) if(color[i] != color[i+1]) {
			for(ui j = i+1;j < vs_size;j ++) assert(color[i] != color[j]);
		}
#endif
		while(true) {
			ui idx = n;
			for(ui i = 0;i < vs_size;i ++) {
				if((i==0||color[i]!=color[i-1])&&(i+1==vs_size||color[i]!=color[i+1])) {
					idx = i;
					break;
				}
			}
			if(idx == n) break;

			current_clique[clique_size ++] = vs[idx];

			ui new_size = 0;
			unsigned char *t_matrix = matrix + vs[idx]*matrix_len;
			for(ui i = 0;i < vs_size;i ++) {
				if(i == idx) continue;
				if(test_bit(t_matrix, vs[i])) {
					vs[new_size] = vs[i];
					color[new_size ++] = color[i];
				}
				else if(i+1 == vs_size||color[i+1] != color[i]) {
					if(new_size == 0||color[i] != color[new_size-1]) return 1;
				}
			}
			vs_size = new_size;
		}

		if(clique_size > max_clique.size()) {
			store_a_larger_clique(clique_size, "kernel_color", 1);
			return 1;
		}
		return 0;
	}

	void move_min_cardinality_color_to_front(ui *vs, const ui vs_size, ui *color) {
		ui c, min_cnt = n, cnt = 0;
		for(ui i = 0;i < vs_size;i ++) {
			++ cnt;
			if(i+1 == vs_size||color[i] != color[i+1]) {
				if(cnt < min_cnt) {
					min_cnt = cnt;
					c = color[i];
				}
				cnt = 0;
			}
		}
		ui idx = vs_size;
		while(idx > 0&&color[idx-1] != c) -- idx;
		for(ui i = idx;i > 0;i --) if(color[i-1] != c) {
			std::swap(vs[i-1], vs[-- idx]);
			std::swap(color[i-1], color[idx]);
		}
	}

	void obtain_degrees(const ui *vs, const ui vs_size, ui *degree) {
		for(ui i = 0;i < vs_size;i ++) degree[vs[i]] = 0;
		for(ui i = 0;i < vs_size;i ++) {
			ui &d = degree[vs[i]];
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
				++ degree[vs[j]];
				++ d;
			}
		}
	}

	char reduce(ui *vs, ui &vs_size, ui *color, const ui threshold, ui clique_size) {
#ifndef NDEBUG
		for(ui i = 0;i < vs_size;i ++) {
			ui d = 0;
			unsigned char *t_matrix = matrix + vs[i]*matrix_len;
			for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) ++ d;
			assert(degree[vs[i]] == d);
		}
#endif

		assert(rid != nullptr);
		ui start = 0;
		for(;start < vs_size;start ++) if(degree[vs[start]] >= threshold) {
			ui *neighbors = rid;
			ui neighbors_n = 0;
			unsigned char *t_matrix = matrix + vs[start]*matrix_len;
			ui color_cnt = 0;
			for(ui i = start+1;i < vs_size;i ++) if(test_bit(t_matrix, vs[i])) {
				neighbors[neighbors_n ++] = i;
				if(!vis[color[i]]) {
					vis[color[i]] = 1;
					++ color_cnt;
				}
			}
			assert(neighbors_n == degree[vs[start]]);
			for(ui i = 0;i < neighbors_n;i ++) {
				vis[color[neighbors[i]]] = 0;
				neighbors[i] = vs[neighbors[i]];
			}
			if(color_cnt < threshold) continue;

			if(degree[vs[start]] >= threshold + 1) break;

			char ok = 1;
			for(ui i = 0;i < neighbors_n&&ok;i ++) {
				t_matrix = matrix + neighbors[i]*matrix_len;
				for(ui j = i+1;j < neighbors_n;j ++) if(!test_bit(t_matrix, neighbors[j])) {
					ok = 0;
					break;
				}
			}
			if(ok) {
				current_clique[clique_size ++] = vs[start];
				for(ui i = 0;i < neighbors_n;i ++) current_clique[clique_size ++] = neighbors[i];
				assert(clique_size == max_clique.size()+1);
				store_a_larger_clique(clique_size, "reduce", 1);
				return 1;
			}
		}

		if(start) {
			for(ui i = start;i < vs_size;i ++) {
				vs[i-start] = vs[i];
				color[i-start] = color[i];
			}
			vs_size -= start;
		}

		if(vs_size == 0) return 1;
		return 0;
	}

	void remove_vertex(ui &current_clique_size, ui u, ui *vs, ui &vs_size, ui *rid, ui *rdegree, std::vector<ui> &del, std::vector<ui> &degree_one, std::vector<ui> &degree_two, std::vector<ui> &degree_three) {
		assert(vs_size > 0&&rid[u] < vs_size);
		vs[rid[u]] = vs[-- vs_size];
		unsigned char *t_matrix = matrix + u*matrix_len;
		for(ui i = 0;i < vs_size;i ++) {
			if(test_bit(t_matrix, vs[i])) {
				if(vs_size - rdegree[vs[i]] + current_clique_size == max_clique.size()) del.pb(vs[i]);
			}
			else {
				-- rdegree[vs[i]];
				switch(rdegree[vs[i]]) {
				case 0: current_clique[current_clique_size ++] = vs[i];
						vs[i] = vs[-- vs_size]; rid[vs[i]] = i;
						-- i;
						break;
				case 1: degree_one.pb(vs[i]);
						break;
				case 2: degree_two.pb(vs[i]);
						break;
				case 3: degree_three.pb(vs[i]);
				}
			}
		}
	}

	void remove_vertex_idx(ui *vs, ui idx, ui &vs_size, ui *rdegree) {
		assert(idx >= 0&&idx < vs_size);

		unsigned char *t_matrix = matrix + vs[idx]*matrix_len;
		vs[idx] = vs[-- vs_size];
		for(ui j = 0;j < vs_size;j ++) if(!test_bit(t_matrix, vs[j])) -- rdegree[vs[j]];
	}

	void search_triangle_matrix(const ui *vs, const ui vs_size, ui &clique_size) {
		if(clique_size+2 == max_clique.size()) {
			for(ui i = 0;i < vs_size&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix1 = matrix + vs[i]*matrix_len;
				assert(rid != nullptr);
				ui *neighbors = rid;
				ui neighbors_n = 0;
				for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix1, vs[j])) {
					neighbors[neighbors_n ++] = vs[j];
				}
				for(ui j = 0;j < neighbors_n&&clique_size < max_clique.size();j ++) {
					unsigned char *t_matrix2 = matrix + neighbors[j]*matrix_len;
					for(ui k = j+1;k < neighbors_n;k ++) if(test_bit(t_matrix2, neighbors[k])) {
						current_clique[clique_size ++] = vs[i];
						current_clique[clique_size ++] = neighbors[j];
						current_clique[clique_size ++] = neighbors[k];
						break;
					}
				}
			}
		}
		else if(clique_size+1 == max_clique.size()) {
			for(ui i = 0;i < vs_size&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix = matrix + vs[i]*matrix_len;
				for(ui j = i+1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
					current_clique[clique_size ++] = vs[i];
					current_clique[clique_size ++] = vs[j];
					break;
				}
			}
		}
		else if(clique_size == max_clique.size()) current_clique[clique_size ++] = vs[0];
	}

	void search_triangle_matrix_color(const ui *vs, const ui vs_size, const ui *color, ui &clique_size) {
		if(clique_size+2 == max_clique.size()) {
			ui idx1 = 1;
			while(color[idx1] == color[0]) ++ idx1;
			assert(idx1+1 < vs_size);
			ui idx2 = idx1+1;
			while(color[idx2] == color[idx1]) ++ idx2;
			assert(idx2 < vs_size);
			for(ui i = 0;i < idx1&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix1 = matrix + vs[i]*matrix_len;
				assert(rid != nullptr);
				ui *neighbors = rid;
				ui neighbors_n = 0;
				for(ui j = idx2;j < vs_size;j ++) if(test_bit(t_matrix1, vs[j])) {
					neighbors[neighbors_n ++] = vs[j];
				}
				for(ui j = idx1;j < idx2&&clique_size < max_clique.size();j ++) if(test_bit(t_matrix1, vs[j])) {
					unsigned char *t_matrix2 = matrix + vs[j]*matrix_len;
					for(ui k = 0;k < neighbors_n;k ++) if(test_bit(t_matrix2, neighbors[k])) {
						current_clique[clique_size ++] = vs[i];
						current_clique[clique_size ++] = vs[j];
						current_clique[clique_size ++] = neighbors[k];
						break;
					}
				}
			}
		}
		else if(clique_size+1 == max_clique.size()) {
			ui idx1 = 1;
			while(color[idx1] == color[0]) ++ idx1;
			assert(idx1 < vs_size);
			for(ui i = 0;i < idx1&&clique_size < max_clique.size();i ++) {
				unsigned char *t_matrix = matrix + vs[i]*matrix_len;
				for(ui j = idx1;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
					current_clique[clique_size ++] = vs[i];
					current_clique[clique_size ++] = vs[j];
					break;
				}
			}
		}
		else if(clique_size == max_clique.size()) current_clique[clique_size ++] = vs[0];
	}

	ui split_vs(ui *vs, const ui vs_size, ui *color, const ui threshold) {
		for(ui i = 0;i < threshold;i ++) head[i] = n;
		for(ui i = vs_size;i > 0;i --) if(color[i-1] < threshold) {
			ui j = i-1;
			next[vs[j]] = head[color[j]];
			head[color[j]] = vs[j];
		}
		ui end_idx = 0;
		for(ui i = 0;i < vs_size;i ++) if(color[i] >= threshold) {
			color[end_idx] = color[i];
			vs[end_idx] = vs[i];
			++ end_idx;
		}
		ui new_size = end_idx;
		for(ui i = threshold;i > 0;i --) for(ui u = head[i-1];u != n;u = next[u]) {
			vs[new_size] = u;
			color[new_size ++] = i-1;
		}
		assert(new_size == vs_size);

		return end_idx;
	}

	ui split_vs_on_color(ui *vs, const ui vs_size, ui *color, const ui threshold) {
		ui max_color = 0;
		for(ui i = 0;i < vs_size;i ++) if(color[i] > max_color) max_color = color[i];
		for(ui i = 0;i <= max_color;i ++) head[i] = n;
		for(ui i = vs_size;i > 0;i --) {
			ui j = i-1;
			next[vs[j]] = head[color[j]];
			head[color[j]] = vs[j];
		}

		ui new_size = 0;
		for(ui i = max_color+1;i > 0;i --) for(ui u = head[i-1];u != n;u = next[u]) {
			vs[new_size] = u;
			color[new_size ++] = i-1;
		}
		assert(new_size == vs_size);

		ui end_idx = 0;
		while(end_idx < vs_size&&color[end_idx] >= threshold) ++ end_idx;

		return end_idx;
	}

	void store_a_larger_clique(const ui clique_size, const char *info, char print) {
		assert(max_clique.size() < clique_size);
		while(max_clique.size() < clique_size) max_clique.pb(0);

		ui contract_size = contractions.size();
		for(ui k = clique_size-1;k > 0;k --) {
			if(current_clique[k] == n) {
				assert(contract_size > 0);
				if(contractions[contract_size-1] == 1) {
					contract_size -= 4;
					if(vis[contractions[contract_size]]) max_clique[k] = contractions[contract_size+1];
					else max_clique[k] = contractions[contract_size+2];
				}
				else if(contractions[contract_size-1] == 2) {
					contract_size -= 5;
					assert(!vis[contractions[contract_size]]||!vis[contractions[contract_size+1]]);
					if(!vis[contractions[contract_size]]&&!vis[contractions[contract_size+1]]) max_clique[k] = contractions[contract_size+3];
					else max_clique[k] = contractions[contract_size+2];
				}
				else if(contractions[contract_size-1] == 3) {
					contract_size -= 5;
					char in1 = (vis[contractions[contract_size]] != 0);
					char in2 = (vis[contractions[contract_size+1]] != 0);
					char in3 = (vis[contractions[contract_size+2]] != 0);
					assert(in1+in2+in3 < 3);
					if(in1+in2+in3 == 0) max_clique[k] = contractions[contract_size+3];
					else if(in1+in2+in3 == 1) {
						if(in1) max_clique[k] = contractions[contract_size+1];
						else if(in2) max_clique[k] = contractions[contract_size+2];
						else max_clique[k] = contractions[contract_size];
					}
					else {
						if(!in1) max_clique[k] = contractions[contract_size];
						else if(!in2) max_clique[k] = contractions[contract_size+1];
						else max_clique[k] = contractions[contract_size+2];
					}
				}
				else {
					printf("WA!\n");
				}
			}
			else max_clique[k] = current_clique[k];
			vis[max_clique[k]] = 1;
		}
		assert(contract_size == 0);
		for(ui k = 1;k < max_clique.size();k ++) vis[max_clique[k]] = 0;
		for(ui k = 1;k < max_clique.size();k ++) max_clique[k] = mapping[max_clique[k]];
		max_clique[0] = current_clique[0];
		if(print) printf("%s finds clique of size: %lu\n", info, max_clique.size());

#ifndef NDEBUG
		ui total_edges = max_clique.size();
		for(ui i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 1;
		for(ui i = 0;i < max_clique.size();i ++) for(ept j = pstart_o[max_clique[i]];j < pstart_o[max_clique[i]+1];j ++) {
			if(vis[edges_o[j]]) total_edges += 2;
		}
		if(total_edges != max_clique.size()*max_clique.size()) printf("WA! Not a clique! %s\n", info);
		for(ui i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 0;
#endif
	}
};

#endif
