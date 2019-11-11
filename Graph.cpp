#include "Utility.h"
#include "Timer.h"
#include "Graph.h"

using namespace std;

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = m = 0;

	pstart = nullptr;
	edges = nullptr;
	pend = nullptr;

	pstart_o = nullptr;
	edges_o = nullptr;

	tmp_vs = nullptr;
	tmp_color = nullptr;

	matrix = nullptr;

	head = nullptr;
	next = nullptr;

	max_core = 0;

	vis = nullptr;

	degree = nullptr;
	rid = nullptr;
	id = nullptr;

	current_clique = nullptr;

	mapping = nullptr;
	mapping_n = matrix_len = 0;

	branches = 0;
}

Graph::~Graph() {
	if(pstart != nullptr) {
		delete[] pstart;
		pstart = nullptr;
	}
	if(edges != nullptr) {
		delete[] edges;
		edges = nullptr;
	}
	if(pend != nullptr) {
		delete[] pend;
		pend = nullptr;
	}
	if(pstart_o != nullptr) {
		delete[] pstart_o;
		pstart_o = nullptr;
	}
	if(edges_o != nullptr) {
		delete[] edges_o;
		edges_o = nullptr;
	}
	if(tmp_vs != nullptr) {
		delete[] tmp_vs;
		tmp_vs = nullptr;
	}
	if(tmp_color != nullptr) {
		delete[] tmp_color;
		tmp_color = nullptr;
	}
	if(matrix != nullptr) {
		delete[] matrix;
		matrix = nullptr;
	}
	if(head != nullptr) {
		delete[] head;
		head = nullptr;
	}
	if(next != nullptr) {
		delete[] next;
		next = nullptr;
	}
	if(vis != nullptr) {
		delete[] vis;
		vis = nullptr;
	}
	if(degree != nullptr) {
		delete[] degree;
		degree = nullptr;
	}
	if(rid != nullptr) {
		delete[] rid;
		rid = nullptr;
	}
	if(id != nullptr) {
		delete[] id;
		id = nullptr;
	}
	if(current_clique != nullptr) {
		delete[] current_clique;
		current_clique = nullptr;
	}
	if(mapping != nullptr) {
		delete[] mapping;
		mapping = nullptr;
	}
}

void Graph::read_graph_binary() {
#ifndef NDEBUG
	printf("\tstart reading graph\n");
#endif
	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("m not equal sum of degrees\n");
#endif

	fclose(f);

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == nullptr) pstart = new ept[n+1];
	if(edges == nullptr) edges = new ui[m];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) {
			fread(edges+pstart[i], sizeof(int), degree[i], f);

			// remove self loops and parallel edges
			ui *buff = edges+pstart[i];
			sort(buff, buff+degree[i]);
			ui idx = 0;
			for(ui j = 0;j < degree[i];j ++) {
				if(buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);
				if(buff[j] == i||(j > 0&&buff[j] == buff[j-1])) continue;
				buff[idx ++] = buff[j];
			}
			degree[i] = idx;
		}

		pstart[i+1] = pstart[i] + degree[i];
	}

	fclose(f);

	delete[] degree;
}

void Graph::read_graph(const char *input_file) {
	printf("# Start reading graph %s\n", input_file);
	FILE *f = Utility::open_file(input_file, "r");

	vector<pair<ui,ui> > vp;
	ui a, b;
	while(fscanf(f, "%u%u", &a, &b) == 2) {
		vp.pb(mp(a,b));
		vp.pb(mp(b,a));
	}
	sort(vp.begin(), vp.end());

	n = 0; m = vp.size();
	for(ui i = 0;i < vp.size();i ++) {
		if(vp[i].first > n) n = vp[i].first;
		if(vp[i].second > n) n = vp[i].second;
		if(i > 0&&vp[i].first == vp[i-1].first&&vp[i].second == vp[i-1].second) printf("WA in read graph\n");
	}
	++ n;

	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	if(pstart != nullptr) delete[] pstart;
	pstart = new ept[n+1];
	if(edges != nullptr) delete[] edges;
	edges = new ui[m];

	pstart[0] = 0;
	ui idx = 0;
	for(ui i = 0;i < n;i ++) {
		pstart[i+1] = pstart[i];
		while(idx < vp.size()&&vp[idx].first == i) edges[pstart[i+1] ++] = vp[idx ++].second;
	}

	fclose(f);

#ifndef NDEBUG
	printf("Finished reading graph\n");
#endif
}

void Graph::heuristic_maximal_clique(char opt) {
	Timer t;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *color = new ui[n];

	max_clique.clear();
	ui UB = degeneracy_maximal_clique_adjacency_list(peel_sequence, core, color, opt);

	delete[] color;
	delete[] core;
	delete[] peel_sequence;

	if(max_clique.size() < UB) printf("\tHeuristic Clique Size: %lu, UB: %u, Total Time: %s (microseconds)\n", max_clique.size(), UB, Utility::integer_to_string(t.elapsed()).c_str());
	else printf("\tMax Clique Size: %lu, Total Time: %s (microseconds)\n", max_clique.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

// exact = 1: compute exact maximum clique
// initial_by_MC_EGO = 1: initialize by MC_EDGE
void Graph::maximum_clique_color(char exact, char initial_by_MC_EGO) {
	Timer t;
	ui *peel_sequence = new ui[n];
	ui *core = new ui[n];
	ui *color = new ui[n];

	max_clique.clear();
	char opt = 1;
	if(!initial_by_MC_EGO) opt = 0;
	ui UB = degeneracy_maximal_clique_adjacency_list(peel_sequence, core, color, opt);
	assert(max_clique.size() >= 2);
	//printf("max_core: %u, UB: %u\n", max_core, UB);
	if(max_clique.size() < UB&&(exact||initial_by_MC_EGO)) {
		ui old_size = max_clique.size();

		ui *out_mapping = nullptr;
		shrink_graph(peel_sequence, core, color, out_mapping, pstart_o, edges_o);
		assert(n > 0);

		if(pend != nullptr) delete[] pend;
		pend = new ept[n];

		assert(current_clique == nullptr); current_clique = new ui[UB];
		assert(mapping == nullptr); mapping = new ui[max_core];
		assert(head == nullptr); head = new ui[max_core];
		assert(next == nullptr); next = new ui[max_core];
		assert(id == nullptr); id = new ui[max_core];
		assert(matrix == nullptr);
#ifdef _BITSET_
		matrix = new unsigned char[max_core*((max_core+7)/8)];
#else
		matrix = new unsigned char[max_core*max_core];
#endif

		assert(degree == nullptr); degree = new ui[n];
		assert(rid == nullptr); rid = new ui[n];
		assert(vis == nullptr); vis = new char[n];
		memset(vis, 0, sizeof(char)*n);

		ui *local_UBs = new ui[n];
		if(initial_by_MC_EGO) {
			UB = ego_degen(peel_sequence, core, color, local_UBs, UB);
			printf("\tMC-EGO Time: %s\n", Utility::integer_to_string(t.elapsed()).c_str());
		}
		else {
			for(ui i = 0;i < n;i ++) local_UBs[i] = n;
		}

		if(exact&&UB > max_clique.size()) search_oriented(peel_sequence, core, color, local_UBs);

		if(max_clique.size() > old_size) {
			for(ui i = 0;i < max_clique.size();i ++) {
				assert(max_clique[i] < n);
				max_clique[i] = out_mapping[max_clique[i]];
			}
		}

		delete[] out_mapping;
		delete[] local_UBs;
	}

	delete[] color;
	delete[] core;
	delete[] peel_sequence;

	if(exact||UB <= max_clique.size()) printf("\tMaximum Clique Size: %lu, Max Depth: %u, Total Time: %s\n", max_clique.size(), max_depth, Utility::integer_to_string(t.elapsed()).c_str());
	else printf("\tHeuristic Clique Size: %lu, UB: %u, Total Time: %s\n", max_clique.size(), UB, Utility::integer_to_string(t.elapsed()).c_str());
}

void Graph::output_clique() {
	FILE *fout = Utility::open_file("clique.txt", "w");
	fprintf(fout, "%lu\n", max_clique.size());
	sort(max_clique.begin(), max_clique.end());
	for(ui i = 0;i < max_clique.size();i ++) fprintf(fout, "%u\n", max_clique[i]);
	fclose(fout);
}

void Graph::verify_clique() {
	vector<ui> clique;
	FILE *fin = Utility::open_file("clique.txt", "r");
	ui clique_n;
	fscanf(fin, "%u", &clique_n);
	for(ui i = 0;i < clique_n;i ++) {
		ui tmp;
		fscanf(fin, "%u", &tmp);
		clique.pb(tmp);
	}
	fclose(fin);

	if(vis == nullptr) vis = new char[n];
	memset(vis, 0, sizeof(char)*n);
	for(ui i = 0;i < clique.size();i ++) {
		if(vis[clique[i]]) {
			printf("WA Clique! Duplicate vertex\n");
			return ;
		}
		vis[clique[i]] = 1;
	}
	for(ui i = 0;i < clique.size();i ++) {
		ui d = 0;
		for(ui j = pstart[clique[i]];j < pstart[clique[i]+1];j ++) if(vis[edges[j]]) ++ d;
		if(d != clique.size()-1) {
			printf("WA Clique! Not adjacent\n");
			return ;
		}
	}
	printf("Correct Clique!\n");
}

ui Graph::coloring_adj_list(const ui *vs, const ui vs_size, const ui original_size, ui *color, char *vis, const ui start_idx, const ui start_color) {
	assert(start_color <= vs_size&&original_size >= vs_size);
	for(ui i = 0;i < original_size;i ++) color[vs[i]] = n;
	for(ui i = vs_size - start_color;i < vs_size;i ++) color[vs[i]] = vs_size - i - 1;

	ui max_color = 0;
	for(ui i = vs_size - start_color;i > start_idx;i --) {
		ui u = vs[i-1];
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui c = color[edges[j]];
			if(c != n) vis[c] = 1;
		}
		for(ui j = 0;;j ++) if(!vis[j]) {
			color[u] = j;
			if(j > max_color) max_color = j;
			break;
		}
		for(ui j = pstart[u];j < pend[u];j ++) {
			ui c = color[edges[j]];
			if(c != n) vis[c] = 0;
		}
	}

	return max_color + 1;
}

ui Graph::degeneracy_ordering_and_coloring_adj(ui *vs, ui *color, const ui vs_size, char *vis, ui *degree, ui *label, const ui K, const ui threshold) {
	for(ui i = 0;i < vs_size;i ++) {
		ui &d = degree[vs[i]];
		d = 0;
		for(ui j = pstart[vs[i]];j < pend[vs[i]]&&label[edges[j]] == K;j ++) ++ d;
	}

	ui start_color = 0;
	for(ui j = 0;j < vs_size;j ++) vis[vs[j]] = 1;
	for(ui j = 0;j < vs_size;j ++) {
		ui min_idx = j;
		for(ui k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
		if(min_idx != j) swap(vs[min_idx], vs[j]);
		ui v = vs[j];
		if(degree[v] + j + 1 == vs_size) {
			start_color = degree[v] + 1;
			break;
		}
		vis[v] = 0;
		for(ui k = pstart[v];k < pend[v]&&label[edges[k]] == K;k ++) if(vis[edges[k]] == 1) -- degree[edges[k]];
	}
	for(ui j = 0;j < vs_size;j ++) vis[vs[j]] = 0;

	for(ui i = 0;i < vs_size;i ++) rid[vs[i]] = i;
	for(ui i = vs_size - start_color;i < vs_size;i ++) color[i] = vs_size - i - 1;

	ui *color_idx = next;

	for(ui i = vs_size - start_color;i > 0;i --) {
		for(ui j = pstart[vs[i-1]];j < pend[vs[i-1]]&&label[edges[j]] == K;j ++) {
			if(rid[edges[j]] >= i) vis[color[rid[edges[j]]]] = 1;
		}
		for(ui j = 0;;j ++) if(!vis[j]) {
			color[i-1] = j;
			break;
		}
		for(ui j = i;j < vs_size;j ++) vis[color[j]] = 0;

		if(color[i-1] < threshold) continue;

		for(ui j = 0;j < threshold;j ++) color_idx[j] = 0;
		for(ui j = pstart[vs[i-1]];j < pend[vs[i-1]]&&label[edges[j]] == K;j ++) if(rid[edges[j]] >= i) {
			ui c = color[rid[edges[j]]];
			if(c < threshold) {
				if(color_idx[c] == 0) color_idx[c] = rid[edges[j]];
				else if(color_idx[c] != n) color_idx[c] = n;
			}
		}

		for(ui j = 0;j < threshold&&color[i-1] >= threshold;j ++) if(color_idx[j] > 0&&color_idx[j] < n) {
			ui v = vs[color_idx[j]];
			for(ui k = pstart[v];k < pend[v]&&label[edges[k]] == K;k ++) {
				if(rid[edges[k]] >= i) vis[color[rid[edges[k]]]] = 1;
			}

			for(ui ii = threshold;ii > 0;) {
				-- ii;
				if(vis[ii]||ii == j) continue;

				color[color_idx[j]] = ii;
				color[i-1] = j;
				break;
			}

			for(ui k = pstart[v];k < pend[v]&&label[edges[k]] == K;k ++) {
				if(rid[edges[k]] >= i) vis[color[rid[edges[k]]]] = 0;
			}
		}
	}

	ui max_color = 0;
	for(ui i = 0;i < vs_size;i ++) if(color[i] > max_color) max_color = color[i];

	return max_color + 1;
}

// degeneracy based clique, if opt also greedy search by maximum degree
// return an upper bound of the maximum clique
ui Graph::degeneracy_maximal_clique_adjacency_list(ui *peel_sequence, ui *core, ui *color, char opt, char greedy_extend) {
	if(opt) heuristic_max_clique_max_degree(10);
	ui threshold = max_clique.size();

	Timer t;

	ui *id_s = peel_sequence;
	ui *degree = new ui[n];
	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	if(pend == nullptr) pend = new ept[n];
	for(ui i = 0;i < n;i ++) {
		pend[i] = pstart[i+1];
		degree[i] = pend[i] - pstart[i];
	}

	ui queue_n = 0, new_size = 0;
	for(ui i = 0;i < n;i ++) if(degree[i] < threshold) id_s[queue_n ++] = i;
	for(ui i = 0;i < queue_n;i ++) {
		ui u = id_s[i]; degree[u] = 0;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(degree[edges[j]] > 0) {
			if((degree[edges[j]] --) == threshold) id_s[queue_n ++] = edges[j];
		}
	}
	for(ui i = 0;i < n;i ++) {
		if(degree[i] >= threshold) id_s[queue_n + (new_size ++)] = i;
		else {
			vis[i] = 1;
			core[i] = 0;
		}
	}
	assert(queue_n + new_size == n);

	ui UB = n;
	if(new_size == 0) UB = max_clique.size();
	else {
		ListLinearHeap *heap = new ListLinearHeap(n, new_size-1);
		heap->init(new_size, new_size-1, id_s+queue_n, degree);
		max_core = 0;
		vector<ui> res;
		for(ui i = 0;i < new_size;i ++) {
			ui u, key;
			heap->pop_min(u, key);
			if(key > max_core) max_core = key;
			core[u] = max_core;
			id_s[queue_n + i] = u;
			if(key + i + 1 == new_size) {
				ui x_size = i+1;
				heap->get_ids(id_s+queue_n, x_size);
				assert(x_size == new_size);
				for(ui j = i;j < new_size;j ++) {
					core[id_s[queue_n+j]] = max_core;
					res.pb(id_s[queue_n+j]);
				}
				break;
			}
			vis[u] = 1;

			for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 0) {
				heap->decrement(edges[j], 1);
			}
		}
		delete heap;
		printf("*** Degeneracy clique size: %lu, max_core: %u, Time: %s (microseconds)\n", res.size(), max_core, Utility::integer_to_string(t.elapsed()).c_str());

		if(res.size() > max_clique.size()) max_clique = res;

		if(max_clique.size() == max_core+1) UB = max_core + 1;
		else{
			memset(vis, 0, sizeof(char)*n);
			ui start_idx = 0;
			while(start_idx < n&&core[id_s[start_idx]] < max_clique.size()) ++ start_idx;
			ui num_color = coloring_adj_list(id_s, n, n, color, vis, start_idx, 0);
			assert(num_color <= UB);
			UB = num_color;
			if(max_clique.size() < UB&&greedy_extend) {
				for(ui i = 0;i < res.size();i ++) vis[res[i]] = 1;
				for(ui i = n-res.size();i > 0;i --) {
					ui u = id_s[i-1], cnt = 0;
					if(core[u] < res.size()) break;
					for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]]) ++ cnt;
					if(cnt == res.size()) {
						res.pb(u);
						vis[u] = 1;
					}
				}
				printf("*** Degen_greedy_extend clique size: %lu, num_colors: %u, Time: %s (microseconds)\n", res.size(), num_color, Utility::integer_to_string(t.elapsed()).c_str());
			}
			if(res.size() > max_clique.size()) max_clique = res;
		}
	}
	memset(vis, 0, sizeof(char)*n);

	delete[] degree;
	delete[] vis;

	return UB;
}

ui Graph::degeneracy_maximal_clique_adj_list(const ui current_clique_size, ui *vs, const ui vs_size, char *vis, ui *degree, ListLinearHeap *heap) {
	assert(vs_size > 0);

	for(ui j = 0;j < vs_size;j ++) vis[vs[j]] = 1;

	char sparse = 0;
	ui total_edges = 0;
	for(ui i = 0;i < vs_size;i ++) total_edges += degree[vs[i]];
	if(total_edges*10 < vs_size*(vs_size-1)) sparse = 1;

	if(sparse) heap->init(vs_size, vs_size-1, vs, degree);

	ui start_color = 0;
	for(ui j = 0;j < vs_size;j ++) {
		ui v, key;
		if(sparse) {
			heap->pop_min(v, key);
			vs[j] = v;
		}
		else {
			ui min_idx = j;
			for(ui k = j+1;k < vs_size;k ++) if(degree[vs[k]] < degree[vs[min_idx]]) min_idx = k;
			if(min_idx != j) swap(vs[min_idx], vs[j]);
			v = vs[j]; key = degree[v];
		}
		if(key + j + 1 == vs_size) {
			start_color = key + 1;
			if(sparse) {
				ui new_size = j+1;
				heap->get_ids(vs, new_size);
				assert(new_size == vs_size);
			}
			if(key + 1 + current_clique_size > max_clique.size()) {
				//printf("Find clique of size %u after search %u egos\n", degree[v] + 2, n-i);
				max_clique.clear();
				max_clique.reserve(key + 1 + current_clique_size);
				for(ui k = j;k < vs_size;k ++) max_clique.pb(vs[k]);
				for(ui k = current_clique_size;k > 0;k --) max_clique.pb(current_clique[k-1]);

#ifndef NDEBUG
				ui total_edges = max_clique.size();
				for(ui i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 1;
				for(ui i = 0;i < max_clique.size();i ++) for(ept j = pstart_o[max_clique[i]];j < pstart_o[max_clique[i]+1];j ++) {
					if(vis[edges_o[j]]) total_edges += 2;
				}
				if(total_edges != max_clique.size()*max_clique.size()) printf("WA! Not a clique\n");
				for(ui i = 0;i < max_clique.size();i ++) vis[max_clique[i]] = 0;
#endif
			}
			break;
		}
		vis[v] = 0;
		if(sparse) {
			for(ui k = pstart[v];k < pend[v];k ++) if(vis[edges[k]] == 1) heap->decrement(edges[k], 1);
		}
		else {
			for(ui k = pstart[v];k < pend[v];k ++) if(vis[edges[k]] == 1) -- degree[edges[k]];
		}
	}
	for(ui j = 0;j < vs_size;j ++) vis[vs[j]] = 0;

	return start_color;
}

//construct a better maximal clique and a tighter upper bound by considering all ego-networks
ui Graph::ego_degen(const ui *peel_sequence, const ui *core, const ui *color, ui *local_UBs, const ui UB) {
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(!vis[i]);
	assert(max_clique.size() >= 2&&max_clique.size() <= n);
#endif

	Timer t;

	ui max_local_UB = 0, initial_size = max_clique.size();
	ui *queue = new ui[max_core];
	ListLinearHeap *heap = new ListLinearHeap(n, max_core);
	heap->init(0, 0, nullptr, nullptr);
	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];
		local_UBs[i-1] = 0;

		if(n-i+1 <= max_clique.size()) continue;
		if(core[u] < max_clique.size()) break;

		//get N^+(u)
		ui vs_size = 0;
		get_higher_neighbors(u, vs_size, vs_buf, color_buf, pstart_o, edges_o);
		assert(vs_size <= max_core);

		//color-based prune
		if(vs_size < max_clique.size()||color_bound(vs_buf.data(), vs_size, color, vis) < max_clique.size()) continue;

		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 1;

		//test max_clique \cup \{u\}
		if(max_clique.size() > initial_size) greedy_extend(u, vis, 0);

		//construct G[N^+(u)] and reduce by k-core
		ui old_size = vs_size, original_size = old_size;
		construct_induced_subgraph(vs_buf.data(), vs_size, vis, degree, pstart_o, edges_o);
		kcore_reduction(vs_buf.data(), vs_size, vis, degree, max_clique.size()-1, queue);
		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
		if(vs_size < old_size&&color_bound(vs_buf.data(), vs_size, color, vis) < max_clique.size()) continue;

		char kernel = 0; // by default set as 0
		ui current_clique_size = 1; current_clique[0] = u;
		if(kernel) {
			//construct matrix
			ui *rdegree = degree;
			construct_matrix(vs_buf.data(), vs_size, mapping, vis, rdegree, pstart_o, edges_o);

			contractions.clear();
			old_size = vs_size;
			degree_one_two_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, contractions);
			if(max_clique.size() >= UB) break;

			//if(vs_size + current_clique_size <= max_clique.size()) continue;
			if(vs_size < old_size&&current_clique_size + color_bound(vs_buf.data(), vs_size, mapping, color, vis) <= max_clique.size()) continue;

			//degeneracy-based maximal clique
			//sort(vs_buf.begin(), vs_buf.begin()+vs_size);
			for(ui j = 0;j < vs_size;j ++) degree[vs_buf[j]] = vs_size - 1 - rdegree[vs_buf[j]];
			ui start_color = degeneracy_maximal_clique_matrix(current_clique_size, vs_buf.data(), vs_size, degree, 1, 0);
			if(max_clique.size() >= UB) break;

			//color-based upper bound
			//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
			local_UBs[i-1] = current_clique_size + coloring_matrix(vs_buf.data(), vs_size, rid, vis, 0, start_color);
		}
		else {
			//degeneracy-based maximal clique
			ui start_color = degeneracy_maximal_clique_adj_list(current_clique_size, vs_buf.data(), vs_size, vis, degree, heap);

			if(max_clique.size() >= UB) break;

			//color-based upper bound
			local_UBs[i-1] = current_clique_size + coloring_adj_list(vs_buf.data(), vs_size, original_size, rid, vis, 0, start_color);
		}

		if(local_UBs[i-1] > max_local_UB) max_local_UB = local_UBs[i-1];
	}

	delete[] queue;
	delete heap;

	if(max_local_UB > UB) max_local_UB = UB;

	ui new_UB = max_clique.size();
	if(max_local_UB > new_UB) new_UB = max_local_UB;

	printf("*** ego_degen clique size: %lu, UB: %u, Time: %s (microseconds)\n", max_clique.size(), new_UB, Utility::integer_to_string(t.elapsed()).c_str());

	return new_UB;
}

void Graph::heuristic_max_clique_max_degree(ui processed_threshold) {
	Timer t;
	assert(max_clique.empty());
	ui *head = new ui[n];
	ui *next = new ui[n];
	ui *degree = new ui[n];

	ui *vis = new ui[n];
	memset(vis, 0, sizeof(ui)*n);

	int max_degree = 0;
	for(ui i = 0;i < n;i ++) head[i] = n;
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1]-pstart[i];
		if(degree[i] > max_degree) max_degree = degree[i];
		next[i] = head[degree[i]];
		head[degree[i]] = i;
	}

	for(ui processed_vertices = 0;max_degree >= max_clique.size()&&processed_vertices < processed_threshold;processed_vertices ++) {
		ui u = n;
		while(max_degree >= max_clique.size()&&u == n) {
			for(ui v = head[max_degree];v != n;) {
				ui tmp = next[v];
				if(degree[v] == max_degree) {
					u = v;
					head[max_degree] = tmp;
					break;
				}
				else if(degree[v] >= max_clique.size()) {
					next[v] = head[degree[v]];
					head[degree[v]] = v;
				}
				v = tmp;
			}
			if(u == n) {
				head[max_degree] = n;
				-- max_degree;
			}
		}
		if(u == n) break;

		vis[u] = 1;
		for(ui k = pstart[u];k < pstart[u+1];k ++) if(!vis[edges[k]]) -- degree[edges[k]];

		vector<ui> vs;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) vs.pb(edges[j]);

		vector<ui> vs_deg(vs.size());
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 2;
		for(ui j = 0;j < vs.size();j ++) {
			ui v = vs[j], d = 0;
			for(ui k = pstart[v];k < pstart[v+1];k ++) {
				if(vis[edges[k]] == 2) ++ d;
			}
			vs_deg[j] = d;
		}
		for(ui j = 0;j < vs.size();j ++) vis[vs[j]] = 0;

		vector<ui> res; res.pb(u);
		ui vs_size = vs.size();
		while(vs_size > 0&&res.size() + vs_size > max_clique.size()) {
			ui idx = 0;
			for(ui j = 1;j < vs_size;j ++) {
				if(vs_deg[j] > vs_deg[idx]) idx = j;
				else if(vs_deg[j] == vs_deg[idx]&&degree[vs[j]] > degree[vs[idx]]) idx = j;
			}
			u = vs[idx];

			ui new_size = 0;
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) vis[edges[j]] = 2;
			for(ui j = 0;j < vs_size;j ++) if(vis[vs[j]]) {
				if(j != new_size) swap(vs[new_size], vs[j]);
				vs_deg[new_size] = vs_deg[j];
				++ new_size;
			}
			for(ui j = pstart[u];j < pstart[u+1];j ++) if(vis[edges[j]] == 2) vis[edges[j]] = 0;

			res.pb(u);
			for(ui k = 0;k < new_size;k ++) vis[vs[k]] = k+2;
			for(ui j = new_size;j < vs_size;j ++) {
				ui v = vs[j];
				for(ui k = pstart[v];k < pstart[v+1];k ++) {
					if(vis[edges[k]] >= 2) -- vs_deg[vis[edges[k]]-2];
				}
			}
			for(ui k = 0;k < new_size;k ++) vis[vs[k]] = 0;

			vs_size = new_size;
		}

		if(res.size() > max_clique.size()) max_clique = res;
	}

	delete[] vis;
	delete[] head;
	delete[] next;
	delete[] degree;

	printf("*** Heuristic clique size: %lu, time: %s (microseconds)\n", max_clique.size(), Utility::integer_to_string(t.elapsed()).c_str());
}

//Build the oriented graph
void Graph::shrink_graph(ui *&peel_sequence, ui *&core, ui *&color, ui *&out_mapping, ept *&pstart_o, ui *&edges_o) {
	ui *rid = new ui[n];
	for(ui i = 0;i < n;i ++) rid[peel_sequence[i]] = i;

	ui clique_size = max_clique.size();
	assert(pstart_o == nullptr);
	pstart_o = new ept[n+1];
	pstart_o[0] = 0;
	ui cnt = 0;
	for(ui i = 0;i < n;i ++) {
		pstart_o[i+1] = pstart_o[i];
		if(core[i] < clique_size) continue;
		++ cnt;
		ept &pos = pstart_o[i+1];
		for(ept j = pstart[i];j < pstart[i+1];j ++) if(rid[edges[j]] > rid[i]) edges[pos ++] = edges[j];
	}

	assert(out_mapping == nullptr);
	out_mapping = new ui[cnt];

	cnt = 0;
	for(ui i = 0;i < n;i ++) if(core[i] >= clique_size) {
		out_mapping[cnt] = i;
		core[cnt] = core[i];
		color[cnt] = color[i];
		rid[i] = cnt ++;
	}
	ui x = n - cnt;
	for(ui i = x;i < n;i ++) peel_sequence[i-x] = rid[peel_sequence[i]];
	assert(edges_o == nullptr);
	edges_o = new ui[pstart_o[n]];
	for(ui i = 0;i < pstart_o[n];i ++) edges_o[i] = rid[edges[i]];
	for(ui i = 0;i < cnt;i ++) pstart_o[i] = pstart_o[out_mapping[i]];
	pstart_o[cnt] = pstart_o[n];
	n = cnt;

	delete[] rid;

	ui *t_peel_sequence = new ui[n];
	memcpy(t_peel_sequence, peel_sequence, sizeof(ui)*n);
	delete[] peel_sequence; peel_sequence = t_peel_sequence;

	ui *t_core = new ui[n];
	memcpy(t_core, core, sizeof(ui)*n);
	delete[] core; core = t_core;

	ui *t_color = new ui[n];
	memcpy(t_color, color, sizeof(ui)*n);
	delete[] color; color = t_color;

	ept *t_pstart = new ept[n+1];
	memcpy(t_pstart, pstart_o, sizeof(ept)*(n+1));
	delete[] pstart_o; pstart_o = t_pstart;

	printf("\tReduced graph size: |V|=%s, |E|=%s (undirected)\n", Utility::integer_to_string(cnt).c_str(), Utility::integer_to_string(pstart_o[n]).c_str());
}

void Graph::recursive_search_clique_color_with_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size, char kernel) {
	assert(clique_size <= max_clique.size()&&vs_buf.size() == color_buf.size()&&vs_size);
	ui *vs = vs_buf.data() + vs_begin;
	ui *color = color_buf.data() + vs_begin;

#ifndef _KERNEL_
	printf("Wrong invocation of recursive_search_clique_color_with_kernelization!\n");
#endif

#ifdef _STATISTIC_
	++ branches;
	if(level > max_depth) max_depth = level;
#endif

#ifndef NDEBUG
	if(!kernel) {
		ui color_cnt = 1;
		for(ui i = 1;i < vs_size;i ++) if(color[i] != color[i-1]) {
			++ color_cnt;
			for(ui j = i+1;j < vs_size;j ++) assert(color[j] != color[i-1]);
		}
		assert(color_cnt + clique_size == max_clique.size()+1);
	}
	for(ui i = 0;i < vs_size;i ++) for(ui j = 1;j < clique_size;j ++) assert(vs[i] != current_clique[j]);
	for(ui i = 0;i < vs_size;i ++) for(ui j = i+1;j < vs_size;j ++) assert(vs[i] != vs[j]);
	for(ui i = 1;i < clique_size;i ++) for(ui j = i+1;j < clique_size;j ++) assert(current_clique[i] == n || current_clique[i] != current_clique[j]);
#endif

	if(clique_size + 3 > max_clique.size()) {
		if(kernel) search_triangle_matrix(vs, vs_size, clique_size);
		else search_triangle_matrix_color(vs, vs_size, color, clique_size);
		if(clique_size > max_clique.size()) store_a_larger_clique(clique_size, "search_triangle", 1);
		return ;
	}

	kernel = 1;

	ui end_idx = 0, old_max_clique_size = max_clique.size();
	if(kernel) {
		obtain_degrees(vs, vs_size, degree);
		if(level) {
			ui *rdegree = degree;
			for(ui i = 0;i < vs_size;i ++) rdegree[vs[i]] = vs_size-1-degree[vs[i]];
			degree_one_two_three_reduction_with_folding_matrix(clique_size, vs, vs_size, rdegree, rid);
			for(ui i = 0;i < vs_size;i ++) degree[vs[i]] = vs_size-1-rdegree[vs[i]];
			if(clique_size > max_clique.size()) store_a_larger_clique(clique_size, "degree_one_two_three_with_folding", 1);

			if(!vs_size||max_clique.size() > old_max_clique_size) return ;
		}

		//ui cc_cnt = compute_connected_components(vs, vs_size);
		//if(cc_cnt > 1) ++ cc_larger_than_one[level];

		//degeneracy-based maximal clique
		ui start_color = degeneracy_maximal_clique_matrix(clique_size, vs, vs_size, degree, 1, 1);

		if(max_clique.size() > old_max_clique_size) return ;

		ui threshold = max_clique.size() - clique_size;
		//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
		ui color_cnt = coloring_matrix_advanced(vs, vs_size, color, start_color, threshold);
		if(color_cnt <= threshold) return ;

		//reduce to (threshold+1)-core
		if(reduce(vs, vs_size, color, threshold, clique_size)) return ;

		//reorganize vs
#ifdef _BRACH_ON_COLOR_
		end_idx = split_vs_on_color(vs, vs_size, color, threshold);
#else
		end_idx = split_vs(vs, vs_size, color, threshold);
#endif

		if(color_cnt <= threshold+1) kernel = 0;
	}
	else {
		move_min_cardinality_color_to_front(vs, vs_size, color);

		end_idx = 1;
		while(end_idx < vs_size&&color[end_idx] == color[0]) ++ end_idx;
	}

	while(vs_buf.size() < vs_begin+vs_size+vs_size) {
		vs_buf.pb(0);
		color_buf.pb(0);
	}
	vs = vs_buf.data() + vs_begin;
	color = color_buf.data() + vs_begin;

	for(ui i = end_idx;i > 0&&max_clique.size() == old_max_clique_size;i --) {
		ui *tvs = vs + vs_size;
		ui *tcolor = color + vs_size;
		ui tvs_end = 0;
		unsigned char *t_matrix = matrix + vs[i-1]*matrix_len;

		ui upper_bound = 0;
		for(ui j = (kernel?i:end_idx);j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
			assert(vs_buf.size() > vs_begin + vs_size + tvs_end);

			if(tvs_end == 0||color[j] != tcolor[tvs_end-1]) ++ upper_bound;

			tvs[tvs_end] = vs[j];
			tcolor[tvs_end ++] = color[j];
		}
		// Note that a vertex u may not connected to a vertex with color i for i < color(u)
		if(upper_bound + 1 + clique_size <= max_clique.size()) continue;

		current_clique[clique_size] = vs[i-1];
		ui new_clique_size = clique_size + 1;
		ui old_size = tvs_end;

		if(!kernel) {
			assert(upper_bound + clique_size == max_clique.size());
			if(kernelization_color(new_clique_size, tvs, tvs_end, tcolor)) continue;
		}

		ui old_contraction_size = contractions.size(), old_changes_size = changes.size();
		recursive_search_clique_color_with_kernelization(level+1, new_clique_size, vs_begin + vs_size, tvs_end, kernel);
		while(contractions.size() > old_contraction_size) contractions.pop_back();
		while(changes.size() > old_changes_size) {
			std::pair<ui, ui> p = changes.back(); changes.pop_back();
			reverse_bit(matrix+p.first*matrix_len, p.second);
			reverse_bit(matrix+p.second*matrix_len, p.first);
		}
		if(!kernel&&old_size + end_idx == vs_size) break;

		vs = vs_buf.data() + vs_begin;
		color = color_buf.data() + vs_begin;
	}
}

void Graph::recursive_search_clique_color_without_kernelization(const ui level, ui clique_size, ui vs_begin, ui vs_size) {
	assert(clique_size <= max_clique.size()&&vs_buf.size() == color_buf.size()&&vs_size);
	ui *vs = vs_buf.data() + vs_begin;
	ui *color = color_buf.data() + vs_begin;

#ifdef _STATISTIC_
	++ branches;
	if(level > max_depth) max_depth = level;
#endif

#ifndef NDEBUG
	//printf("vs_begin: %u\n", vs_begin);
	for(ui i = 0;i < vs_size;i ++) for(ui j = 1;j < clique_size;j ++) assert(vs[i] != current_clique[j]);
	for(ui i = 0;i < vs_size;i ++) for(ui j = i+1;j < vs_size;j ++) assert(vs[i] != vs[j]);
	for(ui i = 1;i < clique_size;i ++) for(ui j = i+1;j < clique_size;j ++) assert(current_clique[i] == n||current_clique[i] != current_clique[j]);
#endif

	if(clique_size == max_clique.size()) {
		current_clique[clique_size ++] = vs[0];
		store_a_larger_clique(clique_size, "search", 0);
		return ;
	}

	ui old_max_clique_size = max_clique.size();
	obtain_degrees(vs, vs_size, degree);

	//degeneracy_ordering
	ui start_color = degeneracy_maximal_clique_matrix(clique_size, vs, vs_size, degree, 1, 0);
	//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);
	if(max_clique.size() > old_max_clique_size) return ;


	//printf("start_color: %u, vs_size: %u\n", start_color, vs_size);

	ui threshold = max_clique.size() - clique_size;
	ui color_cnt = coloring_matrix_advanced(vs, vs_size, color, start_color, threshold);
	if(color_cnt <= threshold) return ;

	//reorganize vs
	ui end_idx = split_vs(vs, vs_size, color, threshold);

	while(vs_buf.size() < vs_begin+vs_size+vs_size) {
		vs_buf.pb(0);
		color_buf.pb(0);
	}
	vs = vs_buf.data() + vs_begin;
	color = color_buf.data() + vs_begin;

	for(ui i = end_idx;i > 0&&max_clique.size() == old_max_clique_size;i --) {
		ui *tvs = vs + vs_size;
		ui *tcolor = color + vs_size;
		ui tvs_end = 0;
		unsigned char *t_matrix = matrix + vs[i-1]*matrix_len;

		ui upper_bound = 0;
		for(ui j = i;j < vs_size;j ++) if(test_bit(t_matrix, vs[j])) {
			assert(vs_buf.size() > vs_begin + vs_size + tvs_end);

			if(tvs_end == 0||color[j] != tcolor[tvs_end-1]) ++ upper_bound;

			tvs[tvs_end] = vs[j];
			tcolor[tvs_end ++] = color[j];
		}
		// Note that a vertex u may not connected to a vertex with color i for i < color(u)
		if(upper_bound + 1 + clique_size <= max_clique.size()) continue;

		current_clique[clique_size] = vs[i-1];
		ui new_clique_size = clique_size + 1;

		recursive_search_clique_color_without_kernelization(level+1, new_clique_size, vs_begin + vs_size, tvs_end);

		vs = vs_buf.data() + vs_begin;
		color = color_buf.data() + vs_begin;
	}
}

void Graph::search_oriented(const ui *peel_sequence, const ui *core, const ui *color, const ui *local_UBs) {
#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) assert(!vis[i]);
	assert(max_clique.size() >= 2&&max_clique.size() <= n);
#endif

#ifdef _STATISTIC_
	ui matrix_cnt = 0;
	double total_density = 0;
	double min_density = 1;
	long total_kernel_effect = 0;

	ui search_ego_cnt = 0;
	double total_ego_density = 0;
	double min_ego_density = 1;

	branches = 0;
	max_depth = 0;
#endif

	ui initial_size = max_clique.size();
	ui *queue = new ui[max_core];

	for(ui i = n;i > 0;i --) {
		ui u = peel_sequence[i-1];

		if(local_UBs[i-1] <= max_clique.size()) continue;
		if(core[u] < max_clique.size()) break;

		ui old_max_clique_size = max_clique.size();

		//get N^+(u)
		ui vs_size = 0;
		get_higher_neighbors(u, vs_size, vs_buf, color_buf, pstart_o, edges_o);
		assert(vs_size <= max_core);

		//color-based prune
		if(vs_size < max_clique.size()||color_bound(vs_buf.data(), vs_size, color, vis) < max_clique.size()) continue;

		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 1;

		//test max_clique \cup \{u\}
		if(max_clique.size() > initial_size&&greedy_extend(u, vis, 1)) {
			for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
			continue;
		}

#ifdef _KERNEL_
		//construct G[N^+(u)] and reduce by k-core
		ui old_size = vs_size;
		construct_induced_subgraph(vs_buf.data(), vs_size, vis, degree, pstart_o, edges_o);
		kcore_reduction(vs_buf.data(), vs_size, vis, degree, max_clique.size()-1, queue);
		for(ui j = 0;j < vs_size;j ++) vis[vs_buf[j]] = 0;
		if(vs_size < old_size&&color_bound(vs_buf.data(), vs_size, color, vis) < max_clique.size()) continue;
#endif

		//construct matrix
		ui *rdegree = degree;
		construct_matrix(vs_buf.data(), vs_size, mapping, vis, rdegree, pstart_o, edges_o);

#ifdef _STATISTIC_
		++ matrix_cnt;
		ui total_edges = 0;
		for(ui j = 0;j < vs_size;j ++) total_edges += vs_size - 1 - rdegree[vs_buf[j]];
		double density = double(total_edges)/(vs_size*(vs_size-1));
		if(density < min_density) min_density = density;
		total_density += density;
#endif

		ui current_clique_size = 1; current_clique[0] = u;
		contractions.clear(); changes.clear();

#ifdef _KERNEL_
		old_size = vs_size;
		//degree_one_two_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, contractions);
		degree_one_two_three_reduction_with_folding_matrix(current_clique_size, vs_buf.data(), vs_size, rdegree, rid);
		if(current_clique_size > max_clique.size()) store_a_larger_clique(current_clique_size, "outside kernelization", 1);
		total_kernel_effect += old_size - vs_size;
		if(max_clique.size() > old_max_clique_size||!vs_size) continue;

#ifdef _STATISTIC_
		++ search_ego_cnt;
		total_edges = 0;
		for(ui j = 0;j < vs_size;j ++) total_edges += vs_size - 1 - rdegree[vs_buf[j]];
		density = double(total_edges)/(vs_size*(vs_size-1));
		if(density < min_ego_density) min_ego_density = density;
		total_ego_density += density;
#endif

		changes.clear();
		recursive_search_clique_color_with_kernelization(0, current_clique_size, 0, vs_size, 1);
		//recursive_search_clique_color_without_kernelization(current_clique_size, 0, vs_size);
#else
		recursive_search_clique_color_without_kernelization(0, current_clique_size, 0, vs_size);
#endif
	}

	delete[] queue;

#ifdef _STATISTIC_
	if(matrix_cnt == 0) printf("No matrix is constructed!\n");
	else {
		printf("Number of matrix constructed: %s\n", Utility::integer_to_string(matrix_cnt).c_str());
		printf("Average density: %.4lf, min density: %.4lf, average kernel_effect: %.4lf\n", total_density/matrix_cnt, min_density, double(total_kernel_effect)/matrix_cnt);

		printf("Number of egos searched: %s, branches: %s\n", Utility::integer_to_string(search_ego_cnt).c_str(), Utility::integer_to_string(branches).c_str());
		if(search_ego_cnt == 0) search_ego_cnt = 1;
		printf("Average ego_density: %.4lf, min ego_density: %.4lf\n", total_ego_density/search_ego_cnt, min_ego_density);
	}
#endif
}
