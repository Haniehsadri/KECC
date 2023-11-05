import org.w3c.dom.ls.LSOutput;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.util.HashMap;
import java.util.Iterator;

import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

public class ECC {

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

    // This set is to implement forced contraction and it contains the edges that has weight at leas k
    //Set<int[]> contractionSet= new HashSet<>();
    Queue<int[]> contractionList = new LinkedList<>();
    boolean check = false;
    //The graph
    HashMap<Integer, Integer>[] mainAdjacencyList;
    int n = 0; //number of vertices

    static int k; // Number of K
    static int t; // Number of Iterations
    long timeforpickrandom = 0l;

    static HashSet<Integer> workSet3 = new HashSet();
    int[][] weights;

    int m;
    int[] translate;
    int[] weight;
    int[] edge;
    int[] pstart;
    // this variables is equal to number of nodes at first, when a node is deleted this variable will decrese by 1
    int finish;

    public ECC(String filename) throws Exception {
        // read in a sequence of pairs of integers (each in the range 0 to N-1),
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        HashMap<Integer, Integer> M = new HashMap();
        List<Integer> node = new ArrayList<Integer>();
        List<Integer> edgesFirst = new ArrayList<Integer>();
        List<Integer> edgesSecond = new ArrayList<Integer>();
        Set<Integer> noderepeat = new HashSet<>();

        Set<String> repeat = new HashSet<>();

        // Get the file extension
        String fileExtension = "";
        int k = filename.lastIndexOf('.');
        if (k > 0) {
            fileExtension = filename.substring(k + 1);
        }

        while ((line = br.readLine()) != null) {

            String[] lineArr = new String[20];
            if ("csv".equalsIgnoreCase(fileExtension)) {
                lineArr = line.split(",");
            } else if ("txt".equalsIgnoreCase(fileExtension)) {
                lineArr = line.split("\\s+");
            }
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
        System.out.println("n" + n);

        finish = n;
        mainAdjacencyList = new HashMap[n];

        m = edgesFirst.size();
        int[] edges1 = edgesFirst.stream().mapToInt(i -> i).toArray();
        int[] edges2 = edgesSecond.stream().mapToInt(i -> i).toArray();
        translate = new int[n];
        edge = new int[m];
        pstart = new int[n + 1];
        Collections.sort(node);
        int size = edges1.length;
        quickSort(edges1, edges2, 0, size - 1);
        for (int i = 0; i < node.size(); i++) {
            M.put(node.get(i), i);
            translate[i] = node.get(i);
        }

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

    public List<Integer> removeVertexWithDegreeLessThanK(int graph_head[],
            int vertex[],
            int next[],
            int reverse[],
            int pre[],
            int degree[],
            List<Integer> inducedsubgraph) {
        HashSet<Integer> workSet3 = new HashSet();
        for (int node : inducedsubgraph) {
            if (degree[node] < k) {
                workSet3.add(node);
            }
        }

        while (!workSet3.isEmpty()) {
            Integer u = workSet3.iterator().next();
            workSet3.remove(u);

            //delete u from its neighbors

            int index = graph_head[u];
            while (index != -1) {
                int temp = next[index];
                int indexofreverse = reverse[index];

                delete_edge(vertex[index], indexofreverse, graph_head, next, pre);

                degree[vertex[index]] = degree[vertex[index]] - 1;

                index = temp;
            }

            int index1 = graph_head[u];
            while (index1 != -1) {
                int temp = next[index1];
                if (degree[vertex[index1]] < k) {
                    workSet3.add(vertex[index1]);

                }
                index1 = temp;
            }
            graph_head[u] = -1;
            inducedsubgraph.remove(u);

        }
        return inducedsubgraph;
    }

    public void decompose() {

        int pend[] = new int[n];
        for (int i = 0; i < n; i++) {
            pend[i] = pstart[i + 1];
        }
        // make copy from pend for construct a partition graph in each itteration
        int Q[] = new int[n];

        int Q_n = 0;
        int computed[] = new int[n];
        int degree[] = new int[n];

        long startTime = System.currentTimeMillis();
        for (int i = 0; i < n; i++) {
            degree[i] = pend[i] - pstart[i];
            if (degree[i] < k) {
                Q[Q_n++] = i;

                computed[i] = 1;
            }
        }

        k_core_prune(k, Q, Q_n, computed, degree, pend);
        int computedForDegree[] = computed.clone();

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
        int counter = 0;
        int vis[] = new int[n];
        List<Set<List<Integer>>> pi = new ArrayList<>();
        int firsti = 0;
        int m = 0;
        while (true) {
            if (computed[m] == 0 || m == computed.length - 1) {
                firsti = m;
                break;
            } else {
                m++;
            }

        }
        int pend3[] = pend.clone();
        int pend5[] = pend.clone();
        int degree2[] = degree.clone();
        construct_pgraph(firsti, Q, vis, computed, pend3, graph_head, vertex, next, reverse, pre);
        pi.add(ContractAndCut(firsti, graph_head, next, pre, reverse, vertex, degree2, degree, pend5, computedForDegree));
        Set<List<Integer>> piI;

        int h = 1;
        int loop = 0;
        while (h < t) {

            System.out.println(h + "itteration starts");

            int components = 0;
            int new_components = 0;

            List<List<List<Integer>>> itteration2 = new ArrayList<>();

            int previous_components = 0;

            pi.add(h, new HashSet());
            piI = pi.get(h);

            for (List<Integer> Gprim : pi.get(h - 1)) {
                components++;
                int[] degree3 = degree.clone();
                int[] computed2 = new int[n];
                Arrays.fill(computed2, 1);
                int b = Gprim.get(0);
                for (int c1 : Gprim) {
                    computed2[c1] = 0;
                }
                for (int u = 0; u < n; u++) {
                    for (int j = pstart[u]; j < pend[u]; j++) {
                        if (computed2[edge[j]] == 1) {
                            if (computedForDegree[edge[j]] == 0) {
                                degree3[u]--;
                            }
                        }
                    }
                }
                int vis2[] = new int[n];

                int pend2[] = pend.clone();

                construct_pgraph(b, Q, vis2, computed2, pend2, graph_head, vertex, next, reverse, pre);
                int pend4[] = pend.clone();
                piI.addAll(ContractAndCut(b, graph_head, next, pre, reverse, vertex, degree3, degree, pend4, computedForDegree));

            }

            new_components = piI.size();
            System.out.println("number of comonents :" + new_components);
            h++;

        }

        long stopTime = System.currentTimeMillis();
        int y = 1;
        int counterc = 0;
        Set<List<Integer>> itteration = pi.get(pi.size() - 1);

        for (List<Integer> component : itteration) {

            for (int s : component) {
                System.out.printf("%d  ", translate[s]);
            }
            System.out.println();
            counterc++;

        }

        System.out.println();

        System.out.println("itteration" + h);
        System.out.println("counter" + counterc);
        System.out.print("Execution time: hh ");
        System.out.println((stopTime - startTime) / 1000.0);
    }

    public void construct_pgraph(int s,
            int Q[],
            int vis[],
            int computed[],
            int pend[],
            int graph_head[],
            int vertex[],
            int next[],
            int reverse[],
            int pre[]) {
        int cnt = 0;

        for (int i = s; i < n; i++) {
            int u = i;

            if (computed[u] == 0) {
                for (int j = pstart[u]; j < pend[u]; ) {
                    int v = edge[j];

                    if (computed[v] == 1) {
                        int temp;
                        temp = edge[j];
                        edge[j] = edge[--pend[u]];
                        edge[pend[u]] = temp;
                    } else {
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

        }
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

    //updat the contractionlist
    public boolean isEdgeInSet(int[] edge) {
        Iterator<int[]> itr = contractionList.iterator();
        while (itr.hasNext()) {
            int[] nextEdge = itr.next();
            if ((edge[0] == nextEdge[0] && edge[1] == nextEdge[1]) || ((edge[0] == nextEdge[1] && edge[1] == nextEdge[0]))) {
                return true;
            }
        }

        return false;
    }

    public void merge(int graph_head[],
            int u,
            int v,
            int vertex[],
            int next[],
            int reverse[],
            int pre[],
            int degree[],
            UF unionFind,
            Set<List<Integer>> allcomponents) {
        int countU = 0;
        int index1 = graph_head[u];

        while (index1 != -1) {
            int tmp = next[index1];
            countU++;
            index1 = tmp;
        }

        int countV = 0;
        int index2 = graph_head[v];

        while (index2 != -1) {
            int tmp = next[index2];
            countV++;

            index2 = tmp;
        }

        if (countU < countV) {
            int temp1 = v;
            v = u;
            u = temp1;
        }

        int index = graph_head[v];

        while (index != -1) {

            int tmp = next[index];

            if (vertex[index] == u) {

                int indexofreverse = reverse[index];

                delete_edge(u, indexofreverse, graph_head, next, pre);

                degree[u] = degree[u] - 1;

            } else {

                assert (vertex[reverse[index]] == v);
                vertex[reverse[index]] = u;

                add_edge(u, graph_head, index, next, pre);

                degree[u] = degree[u] + 1;
            }

            index = tmp;
        }

        //unionFind.union(u, v);
        // check degree of u is less than k
        if (degree[u] < k) {
            workSet3.add(u);
        }

        graph_head[v] = -1;
        //nodeset.remove(v);
        finish--;
        //mainAdjacencyList.remove(v);

        // computed[v] = 1;

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
            // assert(graph_head[u] == indexofE);

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

    public void delete_vertex_recursively(int u,
            int graph_head[],
            int next[],
            int pre[],
            int reverse[],
            int vertex[],
            UF unionFind,
            Set<List<Integer>> allcomponents,
            int degreeSub[],
            int degreOrigin[]) {
        Set<Integer> workset2 = new HashSet<>();
        int index = graph_head[u];
        while (index != -1) {
            int temp = next[index];
            int indexofreverse = reverse[index];

            delete_edge(vertex[index], indexofreverse, graph_head, next, pre);

            degreeSub[vertex[index]] = degreeSub[vertex[index]] - 1;
            if (degreeSub[vertex[index]] < k) {
                workset2.add(vertex[index]);
            }

            index = temp;
        }

        graph_head[u] = -1;
        List<Integer> originalVerticesContractedTou = unionFind.collectComponent(unionFind.find(u));
        allcomponents.add(originalVerticesContractedTou);

        finish--;
    }

    public static int calculateDegree2(HashMap<Integer, HashMap<Integer, Integer>> graph, int u) {
        return graph.get(u).get(-1);
    }

    // This returns a random edge from given graph
    public static int[] pickRandomEdge2(int i, int graph_head[], int next[], int pre[], int reverse[], int vertex[]) {
        int index;

        Random r = new Random();
        int[] edge = new int[2];
        while (true) {
            edge[0] = ThreadLocalRandom.current().nextInt(0, graph_head.length);

            index = graph_head[edge[0]];

            //dont choose -1
            if (index != -1) {

                while (r.nextInt(2) == 1 && next[index] != -1) {
                    index = next[index];
                }
                edge[1] = vertex[index];
                return edge;
            }
        }
    }

    // This is the Algorithm 1 in paper
    public Set<List<Integer>>
    ContractAndCut(int i,
            int graph_head[],
            int next[],
            int pre[],
            int reverse[],
            int vertex[],
            int degreeSub[],
            int degreOrigin[],
            int pend[],
            int basecompute[]) {
        //vaghti ke roye zirgraph sedash mianim bazam n verex darim dar sorati ke kamtar mishe!
        UF unionFind = new UF(n + 1);

        // List<Integer> originalVerticesContractedTou = unionFind.collectComponent(unionFind.find(2));
        // This is a list of all the output componentsA
        Set<List<Integer>> allcomponents = new HashSet<>();

        // Copy the graph G into graph copyG

        while (!isempty(graph_head) || !workSet3.isEmpty()) {

            if (!workSet3.isEmpty()) {

                while (!workSet3.isEmpty()) {

                    int u = workSet3.iterator().next();
                    //System.out.println("u" +u);
                    workSet3.remove(u);
                    List<Integer> originalVerticesContractedTou = unionFind.collectComponent(unionFind.find(u));
                    int[] computed2 = new int[n];
                    int[] degree3 = degreOrigin.clone();
                    Arrays.fill(computed2, 1);
                    int graph_head1[] = new int[n];
                    int vertex1[] = new int[m];
                    int next1[] = new int[m];
                    int reverse1[] = new int[m];
                    int pre1[] = new int[m];
                    Arrays.fill(graph_head1, -1);
                    Arrays.fill(vertex1, -1);
                    Arrays.fill(next1, -1);
                    Arrays.fill(reverse1, -1);
                    Arrays.fill(pre1, -1);
                    int Q[] = new int[n];

                    //avalin khone component ro migigre
                    int b = originalVerticesContractedTou.get(0);
                    for (int c1 : originalVerticesContractedTou) {
                        computed2[c1] = 0;
                    }

                    for (int y = 0; y < n; y++) {
                        for (int j = pstart[y]; j < pend[y]; j++) {
                            if (computed2[edge[j]] == 1) {
                                if (basecompute[edge[j]] != 1) {
                                    degree3[y] = degree3[y] - 1;
                                }
                            }
                        }
                    }
                    int vis2[] = new int[n];
                    construct_pgraph(b, Q, vis2, computed2, pend, graph_head1, vertex1, next1, reverse1, pre1);
                    originalVerticesContractedTou = removeVertexWithDegreeLessThanK(graph_head1, vertex1, next1, reverse1, pre1, degree3,
                            originalVerticesContractedTou);

                    if (!originalVerticesContractedTou.isEmpty()) {
                        allcomponents.add(originalVerticesContractedTou);
                    }

                    //delete u from its neighbors

                    int index = graph_head[u];
                    //System.out.println("u"+u);
                    while (index != -1) {
                        int temp = next[index];
                        int indexofreverse = reverse[index];

                        delete_edge(vertex[index], indexofreverse, graph_head, next, pre);

                        degreeSub[vertex[index]] = degreeSub[vertex[index]] - 1;

                        index = temp;
                    }

                    int index1 = graph_head[u];
                    while (index1 != -1) {

                        int temp = next[index1];
                        if (degreeSub[vertex[index1]] < k) {
                            workSet3.add(vertex[index1]);
                            if (vertex[index1] == 262010) {
                                System.out.println("line 789"
                                );
                            }
                        }
                        index1 = temp;
                    }

                    //remove u
                    if (u == 262010) {
                        System.out.println("801");
                    }
                    graph_head[u] = -1;

                }
            } else {

                long startTime = System.currentTimeMillis();
                int[] randomEdge = pickRandomEdge2(i, graph_head, next, pre, reverse, vertex);
                long stopTime = System.currentTimeMillis();
                timeforpickrandom += (stopTime - startTime);

                merge(graph_head, randomEdge[0], randomEdge[1], vertex, next, reverse, pre, degreeSub, unionFind, allcomponents);

                unionFind.union(randomEdge[0], randomEdge[1]);
            }
        }
        return allcomponents;
    }

    public boolean isempty(int[] graph_head) {
        for (int i = 0; i < n; i++) {
            if (graph_head[i] != -1) {
                return false;
            }
        }
        return true;
    }

    public static void main(String[] argv) throws Exception {

        // Check for enough command-line arguments
        if (argv.length != 3) {
            System.out.println("Usage: <path_to_file> <k_value> <t_value>");
            return;
        }

        // Parse command-line arguments
        String filePath = argv[0];

        k = Integer.parseInt(argv[1]);
        t = Integer.parseInt(argv[2]);

        // Initialize ECC with given file path
        ECC ecc = new ECC(filePath);

        ecc.decompose();
    }
}


