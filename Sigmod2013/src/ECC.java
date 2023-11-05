import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.ThreadLocalRandom;

public class ECC {

    static class Pair {
        int u, v;

        public Pair(int u, int v) {
            this.u = u;
            this.v = v;
        }
    }

    static class Sorting implements Comparator<Pair> {
        public int compare(Pair p1, Pair p2) {
            if (p2.u == p1.u) {
                return p1.v - p2.v;
            }
            return p1.u - p2.u;
        }
    }

    static void swap(int[] arr, int i, int j) {
        int temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }

    static int partition(int[] arr1, int[] arr2, int low, int high) {

        // pivot
        int random = ThreadLocalRandom.current().nextInt(low, high + 1);

        swap(arr1, random, high);
        swap(arr2, random, high);
        int pivot = arr1[high];
        int pivot2 = arr2[high];

        // Index of smaller element and
        // indicates the right position
        // of pivot found so far
        int i = (low - 1);

        for (int j = low; j <= high - 1; j++) {

            // If current element is smaller
            // than the pivot
            if (arr1[j] < pivot) {

                // Increment index of
                // smaller element
                i++;
                swap(arr1, i, j);
                swap(arr2, i, j);
            } else if (arr1[j] == pivot && arr2[j] < pivot2) {
                // Increment index of
                // smaller element
                i++;
                swap(arr1, i, j);
                swap(arr2, i, j);
            }
        }
        swap(arr1, i + 1, high);
        swap(arr2, i + 1, high);
        return (i + 1);
    }

    static void quickSort(int[] arr1, int[] arr2, int low, int high) {
        if (low < high) {

            // pi is partitioning index, arr[p]
            // is now at right place
            int pi = partition(arr1, arr2, low, high);

            // Separately sort elements before
            // partition and after partition
            quickSort(arr1, arr2, low, pi - 1);
            quickSort(arr1, arr2, pi + 1, high);
        }
    }

    int n;
    int m;
    int[] translate;
    int[] edge;
    int[] pstart;

    public ECC(String filename) throws Exception {
        // read in a sequence of pairs of integers (each in the range 0 to N-1),
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        HashMap<Integer, Integer> M = new HashMap();
        List<Integer> node = new ArrayList<Integer>();
        List<Integer> edgesFirst = new ArrayList<Integer>();
        List<Integer> edgesSecond = new ArrayList<Integer>();

        Set<String> repeat = new HashSet<>();
        Set<Integer> noderepeat = new HashSet<>();

        // Get the file extension
        String fileExtension = "";
        int k = filename.lastIndexOf('.');
        if (k > 0) {
            fileExtension = filename.substring(k + 1);
        }

        // Choose the delimiter based on the file extension
        while ((line = br.readLine()) != null) {

            String[] lineArr = new String[20];
            if ("csv".equalsIgnoreCase(fileExtension)) {
                lineArr = line.split(",");
            } else if ("txt".equalsIgnoreCase(fileExtension)) {
                lineArr = line.split("\\s+");
            }
            //csv
            int p = Integer.parseInt(lineArr[0]);
            int q = Integer.parseInt(lineArr[1]);
            if (!noderepeat.contains(p)) {
                node.add(p);
                noderepeat.add(p);
            }
            if (!noderepeat.contains(q)) {
                node.add(q);
                noderepeat.add(q);
            }
            if (!repeat.contains(p + "," + q)) {
                edgesFirst.add(p);
                edgesSecond.add(q);

                edgesFirst.add(q);
                edgesSecond.add(p);
                repeat.add(p + "," + q);
                repeat.add(q + "," + p);
            }
        }

        br.close();
        System.out.println("read");

        n = node.size();

        m = edgesFirst.size();
        int[] edges1 = edgesFirst.stream().mapToInt(i -> i).toArray();
        int[] edges2 = edgesSecond.stream().mapToInt(i -> i).toArray();
        translate = new int[n];
        edge = new int[m];
        pstart = new int[n + 1];
        Collections.sort(node);
        int size = edges1.length;
        quickSort(edges1, edges2, 0, size - 1);
        System.out.println("sorted");
        for (int i = 0; i < node.size(); i++) {
            M.put(node.get(i), i);
            translate[i] = node.get(i);
        }
        ;

        char preserved = 1;
        for (int i = 0; i < node.size(); i++) {
            if (node.get(i) != i) {
                preserved = 0;
            }
        }
        if (preserved != 1) {
            System.out.println("Node ids are not preserved!\n");
        }

        int j = 0;
        for (int i = 0; i < n; i++) {
            pstart[i] = j;
            while (j < m && edges1[j] == node.get(i)) {
                edge[j] = M.get(edges2[j]);
                ++j;
            }
        }
        pstart[n] = j;
    }

    public void k_edge_connected_component(int K, boolean print) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.currentTimeMillis();

        if (K < 2) {
            System.out.printf("K must be at least 2!\n");
            return;
        }

        int pend[] = new int[n];
        for (int i = 0; i < n; i++) {
            pend[i] = pstart[i + 1];
        }

        int computed[] = new int[n];
        Arrays.fill(computed, 0);

        int degree[] = new int[n];

        int Q[] = new int[n];
        int Q_n = 0;

        // k-core-based pruning
        for (int i = 0; i < n; i++) {
            degree[i] = pend[i] - pstart[i];
            if (degree[i] < K) {
                Q[Q_n++] = i;

                computed[i] = 1;
            }
        }

        k_core_prune(K, Q, Q_n, computed, degree, pend);

        int ids[] = new int[n];
        int cstart[] = new int[n + 1];

        int c_n = 0;
        int counter = 0;

        cstart[0] = 0;

        int graph_head[] = new int[n];
        int vertex[] = new int[m];
        int next[] = new int[m];
        int reverse[] = new int[m];
        int pre[] = new int[m];

        for (int i = 0; i < m; i++) {
            next[i] = -1;
        }

        for (int i = 0; i < m; i++) {
            reverse[i] = -1;
        }

        for (int i = 0; i < m; i++) {
            pre[i] = -1;
        }

        for (int i = 0; i < n; i++) {
            graph_head[i] = -1;
        }

        int sv_next[] = new int[n];
        int sv_last[] = new int[n];

        int vis[] = new int[n];
        Arrays.fill(vis, 0);

        ListLinearHeap heap = new ListLinearHeap(n, K);

        int keys[] = new int[n];
        int countj = 0;
        for (int i = 0; i < n; ) {
            if (computed[i] == 1) {
                ++i;
                continue;
            }

            boolean dide = false;
            int Qn1 = construct_pgraph(i, Q, vis, computed, pend, graph_head, vertex, next, reverse, pre, sv_next, sv_last);

            int new_cn = decomposition(i, K, cstart, ids, c_n, Q, keys, vis, graph_head, vertex, next, reverse, pre, heap, sv_next,
                    sv_last);

            if (new_cn == c_n + 1) {
                int start = cstart[c_n];
                int end = cstart[c_n + 1];
                for (int j = start; j < end; j++) {
                    computed[ids[j]] = 1;
                }
                ++c_n;
            } else {
                remove_inter_edges(K, c_n, new_cn, cstart, ids, keys, pend, Q, computed, degree);

            }
            countj++;
        }
        long stopTime = System.currentTimeMillis();
        if (print) {
            print_kecc(K, c_n, cstart, ids);
        }
        System.out.print("Execution time: ");
        System.out.println((stopTime - startTime) / 1000.0);
    }

    void k_core_prune(int K, int Q[], int Q_n, int computed[], int degree[], int pend[]) {
        for (int i = 0; i < Q_n; i++) {
            int u = Q[i];
            for (int j = pstart[u]; j < pend[u]; j++) {
                if (computed[edge[j]] == 0) {
                    int v = edge[j];
                    if (degree[v] == K) {
                        Q[Q_n++] = v;
                        computed[v] = 1;
                    }
                    --degree[v];
                }
            }
        }
    }

    public int construct_pgraph(int s,
            int Q[],
            int vis[],
            int computed[],
            int pend[],
            int graph_head[],
            int vertex[],
            int next[],
            int reverse[],
            int pre[],
            int sv_next[],
            int sv_last[]) {
        int cnt = 0;
        Q[0] = s;
        vis[s] = 1;
        int Q_n1 = 1;

        for (int i = 0; i < Q_n1; i++) {
            int u = Q[i];
            sv_next[u] = sv_last[u] = u;
            for (int j = pstart[u]; j < pend[u]; ) {
                int v = edge[j];
                if (computed[v] == 1) {
                    int temp;
                    temp = edge[j];
                    edge[j] = edge[--pend[u]];
                    edge[pend[u]] = temp;
                } else {
                    if (vis[v] == 0) {
                        Q[Q_n1++] = v;
                        vis[v] = 1;
                    }

                    if (v > u) {
                        vertex[cnt] = v;
                        reverse[cnt] = cnt + 1;
                        //add parameter
                        add_edge(u, graph_head, cnt, next, pre);

                        ++cnt;

                        vertex[cnt] = u;
                        reverse[cnt] = cnt - 1;
                        add_edge(v, graph_head, cnt, next, pre);

                        ++cnt;
                    }
                    ++j;
                }
            }
        }
        for (int i = 0; i < Q_n1; i++) {
            vis[Q[i]] = 0;
        }
        return Q_n1;
    }

    public int decomposition(int s,
            int K,
            int cstart[],
            int ids[],
            int c_n,
            int Q[],
            int keys[],
            int vis[],
            int graph_head[],
            int vertex[],
            int next[],
            int reverse[],
            int pre[],
            ListLinearHeap heap,
            int sv_next[],
            int sv_last[]) {
        while (graph_head[s] != -1) {
            heap.init(0, K, null, null);
            heap.insert(s, 0);
            WrapInt u = new WrapInt();
            WrapInt key = new WrapInt();
            u.value = 0;
            key.value = 0;
            int L_n = 0;
            int L[] = new int[n];
            while (heap.pop_max(u, key)) {
                //Q[Q_n++] = u;
                L[L_n++] = u.value;
                vis[u.value] = 1; //vis[u] = 1 means u is in L, vis[u] = 2 means u is in heap, and vis[u] = 3 means u is in Q
                keys[u.value] = key.value;

                //ui new_Qn = Q_n;
                int new_Qn = 0;
                //initilize a q with u

                Q[new_Qn++] = u.value;

                for (int j = 0; j < new_Qn; j++) {

                    int v = Q[j];

                    int index = graph_head[v];

                    while (index != -1) {
                        if (vis[vertex[index]] != 1) {
                            int w = vertex[index];
                            if (vis[w] == 3) {
                                ++keys[w];
                                index = next[index];
                                continue;
                            }

                            if (vis[w] == 2) {
                                key.value = heap.remove(w);
                            } else {
                                key.value = 0;
                            }

                            ++key.value;
                            if (key.value >= K) {
                                Q[new_Qn++] = w;
                                keys[w] = key.value;
                                vis[w] = 3;
                            } else {
                                heap.insert(w, key.value);
                                vis[w] = 2;
                            }
                        }
                        index = next[index];
                    }

                    if (v == u.value) {
                        continue;
                    }

                    // contract u and v
                    vis[v] = 0;
                    keys[u.value] += keys[v];
                    merge(graph_head, u.value, v, keys, sv_next, sv_last, vertex, next, reverse, pre);
                }
            }

            while (L_n > 0 && keys[L[L_n - 1]] < K) {
                int u1 = L[--L_n];

                vis[u1] = 0;

                int index = graph_head[u1];

                while (index != -1) {

                    int indexofreverse = reverse[index];

                    delete_edge(vertex[index], indexofreverse, graph_head, next, pre);
                    index = next[index];
                }
                graph_head[u1] = -1;

                int pos = cstart[c_n];
                ids[pos++] = u1;

                while (sv_next[u1] != u1) {
                    u1 = sv_next[u1];
                    ids[pos++] = u1;
                }
                cstart[++c_n] = pos;

            }

            for (int i = 0; i < L_n; i++) {
                vis[L[i]] = 0;
            }
        }

        return c_n;
    }

    public void merge(int graph_head[],
            int u,
            int v,
            int keys[],
            int sv_next[],
            int sv_last[],
            int vertex[],
            int next[],
            int reverse[],
            int pre[]) {
        int index = graph_head[v];

        while (index != -1) {

            int tmp = next[index];

            if (vertex[index] == u) {
                --keys[u];

                int indexofreverse = reverse[index];

                delete_edge(u, indexofreverse, graph_head, next, pre);
            } else {
                vertex[reverse[index]] = u;

                add_edge(u, graph_head, index, next, pre);
            }

            index = tmp;
        }
        graph_head[v] = -1;

        sv_next[sv_last[u]] = v;
        sv_last[u] = sv_last[v];
    }

    public void add_edge(int u, int graph_head[], int index, int next[], int pre[]) {
        if (graph_head[u] != -1) {
            pre[graph_head[u]] = index;
        }
        next[index] = graph_head[u];
        graph_head[u] = index;
        pre[index] = -1;
    }

    public void delete_edge(int u, int indexofE, int graph_head[], int next[], int pre[]) {

        if (pre[indexofE] == -1) {
            if (next[indexofE] != -1) {
                indexofE = next[indexofE];

                pre[indexofE] = -1;
                graph_head[u] = indexofE;
            } else {
                graph_head[u] = -1;
            }
        } else {
            next[pre[indexofE]] = next[indexofE];
            if (next[indexofE] != -1) {
                pre[next[indexofE]] = pre[indexofE];
            }
        }
    }

    public void remove_inter_edges(int K,
            int c_n,
            int new_cn,
            int cstart[],
            int ids[],
            int component_id[],
            int pend[],
            int Q[],
            int computed[],
            int degree[]) {
        for (int i = c_n; i < new_cn; i++) {
            for (int j = cstart[i]; j < cstart[i + 1]; j++) {
                component_id[ids[j]] = i;
            }
        }

        for (int i = c_n; i < new_cn; i++) {
            for (int j = cstart[i]; j < cstart[i + 1]; j++) {
                int u = ids[j];
                for (int k = pstart[u]; k < pend[u]; ) {
                    int v = edge[k];
                    if (component_id[v] != component_id[u]) {
                        int temp;
                        temp = edge[k];
                        edge[k] = edge[--pend[u]];
                        edge[pend[u]] = temp;

                    }
                    //swap(edges[k], edges[--pend[u]]);
                    else {
                        ++k;
                    }

                }
            }
        }

        int Q_n = 0;
        for (int i = c_n; i < new_cn; i++) {
            for (int j = cstart[i]; j < cstart[i + 1]; j++) {
                int u = ids[j];
                degree[u] = pend[u] - pstart[u];
                if (degree[u] < K) {
                    Q[Q_n++] = u;
                    computed[u] = 1;
                }
            }
        }

        k_core_prune(K, Q, Q_n, computed, degree, pend);
    }

    public void print_kecc(int K, int c_n, int cstart[], int ids[]) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer2 = new PrintWriter("the-file-name5.txt", "UTF-8");
        PrintStream ps = new PrintStream("file_location.txt");
        System.out.println("cn " + c_n);
        for (int i = 0; i < c_n; i++) {
            int start = cstart[i];
            int end = cstart[i + 1];
            Arrays.sort(ids, start, end);
            ;
            for (int j = start; j < end; j++) {
                writer2.printf("%d  ", translate[ids[j]]);
                int k = translate[ids[j]];
                ps.print(k);
                ps.print(" ");
                writer2.printf(" ");
                System.out.printf("%d  ", translate[ids[j]]);
            }
            System.out.println();
            writer2.println();
            ps.println();
        }
        System.out.println("number of components " + c_n);
        writer2.close();
    }

    public static void main(String[] argv) throws Exception {
        if (argv.length != 2) {
            System.out.println("Please provide the file path and the value of k.");
            System.out.println("Usage: java YourClassName <file_path> <k_value>");
            return;
        }

        // Retrieve file path and k value from argv
        String filePath = argv[0];
        int k;

        k = Integer.parseInt(argv[1]);

        // }

        // Create ECC object using the provided file path
        ECC ecc = new ECC(filePath);

        ecc.k_edge_connected_component(6, true);
    }
}
