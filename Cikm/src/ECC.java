import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.Iterator;

import java.util.Set;

public class ECC {
    // This set is to implement forced contraction and it contains the edges that has weight at leas k
    static Queue<int[]> contractionList = new LinkedList();

    //The graph
    static HashMap<Integer, HashMap<Integer, Integer>> mainAdjacencyList = new HashMap();
    static int n = 0; //number of vertices
    static int k = 7; // Number of K
    static int t = 9;
    static PrintWriter writer;
    static StringBuilder sb;

    static HashSet<Integer> workSet = new HashSet();

    public ECC(String filename) throws Exception {
        // read in a sequence of pairs of integers (each in the range 0 to N-1),
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line;
        HashMap<Integer, Integer> neighborsOfP;
        HashMap<Integer, Integer> neighborsOfQ;
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

            if (p > n) {
                n = p;
            }
            if (q > n) {
                n = q;
            }

            //Add edge {p,q}. This translates to:
            //adding (q,1) to adjacency list of p, and
            //adding (p,1) to adjacency list of q.
            if (!mainAdjacencyList.containsKey(p)) {
                mainAdjacencyList.put(p, new HashMap());
                mainAdjacencyList.get(p).put(-1, 0); //Initially the degree of vertex is 0
            }
            // In the bigger file there are pq and qp edges
            if (!mainAdjacencyList.get(p).containsKey(q)) {
                mainAdjacencyList.get(p).put(q, 1); //Initially the edge is weighted 1
                mainAdjacencyList.get(p).put(-1, mainAdjacencyList.get(p).get(-1) + 1); // Add one to degree of vertex after adding an edge
            }

            if (!mainAdjacencyList.containsKey(q)) {
                mainAdjacencyList.put(q, new HashMap());
                mainAdjacencyList.get(q).put(-1, 0); //Initially the degree of vertex is 0
            }
            if (!mainAdjacencyList.get(q).containsKey(p)) {
                mainAdjacencyList.get(q).put(p, 1); //Initially the edge is weighted 1
                mainAdjacencyList.get(q).put(-1, mainAdjacencyList.get(q).get(-1) + 1); // Add one to degree of vertex after adding an edge
            }

        }
        br.close();
        System.out.println("read");

        // add vertices with degree less than k to workset
        Iterator<Integer> iterVertices = mainAdjacencyList.keySet().iterator();
        while (iterVertices.hasNext()) {
            int u = iterVertices.next();
            if (calculateDegree(mainAdjacencyList, u) < k) {
                workSet.add(u);

            }
        }
        //continue till workset gets empty
        while (!workSet.isEmpty()) {
            int u = workSet.iterator().next();
            workSet.remove(u);

            //let's delete u from all it's neighbors adjacency lists
            Object[] neighbors = mainAdjacencyList.get(u).keySet().toArray();
            for (int i = 0; i < neighbors.length; i++) {
                int v = (int) neighbors[i]; //neighbor of u
                if (v == -1) {
                    continue;
                }
                if (v == u) //self-loop
                {
                    continue;
                }
                //delete u from adjacency list of v and reduce the degree of v
                HashMap<Integer, Integer> hashMapOfV = mainAdjacencyList.get(v);
                int weight = hashMapOfV.get(u);
                hashMapOfV.remove(u);
                //set the degree of v
                hashMapOfV.put(-1, hashMapOfV.get(-1) - weight);
                //if degree of v is smaller than k, add v to queue
                if (calculateDegree(mainAdjacencyList, v) < k) {
                    workSet.add(v);
                }
            }

            //now let's remove u from G
            mainAdjacencyList.remove(u);
        }
        System.out.println("hi");
    }

    // This return true if there is a vertex with degree less than k
    public static boolean isThereAnyVertexWithDegreeLessThanK(HashMap<Integer, HashMap<Integer, Integer>> graph) {
        for (Map.Entry<Integer, HashMap<Integer, Integer>> entry : graph.entrySet()) {
            int vertex = entry.getKey();
            if (calculateDegree(graph, vertex) < k) {
                // return true if at least one of the vertices has degree less than k
                return true;
            }
        }
        // return false if none of the vertices has degree less than k
        return false;
    }

    public static HashMap<Integer, HashMap<Integer, Integer>> recursivelyRemoveVertexWithDegreeLessK(HashMap<Integer, HashMap<Integer, Integer>> mainAdjacencyList,
            int K) {
        HashSet<Integer> workSet3 = new HashSet();
        Iterator<Integer> iterVertices = mainAdjacencyList.keySet().iterator();
        while (iterVertices.hasNext()) {
            int u = iterVertices.next();
            if (calculateDegree(mainAdjacencyList, u) < K) {
                workSet3.add(u);

            }
        }
        //continue till workset gets empty
        while (!workSet3.isEmpty()) {
            int u = workSet3.iterator().next();
            workSet3.remove(u);

            //let's delete u from all it's neighbors adjacency lists
            Object[] neighbors = mainAdjacencyList.get(u).keySet().toArray();
            for (int i = 0; i < neighbors.length; i++) {
                int v = (int) neighbors[i]; //neighbor of u
                if (v == -1) {
                    continue;
                }
                if (v == u) //self-loop
                {
                    continue;
                }
                //delete u from adjacency list of v and reduce the degree of v
                HashMap<Integer, Integer> hashMapOfV = mainAdjacencyList.get(v);
                int weight = hashMapOfV.get(u);
                hashMapOfV.remove(u);
                //set the degree of v
                hashMapOfV.put(-1, hashMapOfV.get(-1) - weight);
                //if degree of v is smaller than k, add v to queue
                if (calculateDegree(mainAdjacencyList, v) < k) {
                    workSet3.add(v);
                }
            }

            //now let's remove u from G
            mainAdjacencyList.remove(u);
        }

        return mainAdjacencyList;
    }

    // This calculates the degree of a given vertex
    public static int calculateDegree(HashMap<Integer, HashMap<Integer, Integer>> graph, int u) {
        if (graph.get(u) == null) {
            System.out.println("a");
        }
        return graph.get(u).get(-1);
    }

    // This builds the Induced SubGraph
    public static HashMap<Integer, HashMap<Integer, Integer>> buildInducedSubGraph(Set<Integer> inducedVertices) {
        HashMap<Integer, HashMap<Integer, Integer>> inducedAdjacencyList = new HashMap();
        for (int u : inducedVertices) {
            inducedAdjacencyList.put(u, new HashMap());
            inducedAdjacencyList.get(u).put(-1, 0);

            Set<Integer> neighbours = mainAdjacencyList.get(u).keySet();
            int weight;
            for (int v : neighbours) {
                if (inducedVertices.contains(v)) {
                    //Add edge to AdjacencyList
                    weight = mainAdjacencyList.get(u).get(v);
                    inducedAdjacencyList.get(u).put(v, weight);
                    inducedAdjacencyList.get(u).put(-1, inducedAdjacencyList.get(u).get(-1) + weight);
                }
            }
        }
        return inducedAdjacencyList;
    }

    // This outputs the graph
    public static void output(HashMap<Integer, HashMap<Integer, Integer>> graph) {

        for (Map.Entry<Integer, HashMap<Integer, Integer>> entry : graph.entrySet()) {
            int vertex = entry.getKey();
            System.out.print("  " + vertex + "  ");
        }
        System.out.println();

    }

    public static HashMap<Integer, HashMap<Integer, Integer>> graphDeepCopy(HashMap<Integer, HashMap<Integer, Integer>> adjacencyListG) {

        HashMap<Integer, HashMap<Integer, Integer>> copyAdjacencyList = new HashMap();

        //Copy Adjacency list of graph into a new Adjacency list
        for (Map.Entry<Integer, HashMap<Integer, Integer>> entry : adjacencyListG.entrySet()) {
            int vertex = entry.getKey();
            HashMap<Integer, Integer> value = entry.getValue();
            copyAdjacencyList.put(vertex, new HashMap());
            for (Map.Entry<Integer, Integer> entry2 : value.entrySet()) {
                int neighbour = entry2.getKey();
                int weight = entry2.getValue();
                //Add edg to AdjacencyList of copied Graph
                copyAdjacencyList.get(vertex).put(neighbour, weight);
            }
        }

        return copyAdjacencyList;
    }

    public static boolean isEdgeInSet(int[] edge) {
        Iterator<int[]> itr = contractionList.iterator();
        while (itr.hasNext()) {
            int[] nextEdge = itr.next();
            if ((edge[0] == nextEdge[0] && edge[1] == nextEdge[1]) || ((edge[0] == nextEdge[1] && edge[1] == nextEdge[0]))) {
                return true;
            }
        }

        return false;
    }

    public static void forcedContract(HashMap<Integer, HashMap<Integer, Integer>> graphAdjacencyList, int[] edge) {
        HashMap<Integer, Integer> hashStart = graphAdjacencyList.get(edge[0]);
        HashMap<Integer, Integer> hashEnd = graphAdjacencyList.get(edge[1]);

        int contractedFrom;
        int contractedTo;

        HashMap<Integer, Integer> hashMapSmaller;
        HashMap<Integer, Integer> hashMapBigger;

        // Here we select the smaller hash map that need to be merged into the bigger hash map
        if (hashStart.size() > hashEnd.size()) {
            hashMapSmaller = hashEnd;
            hashMapBigger = hashStart;
            contractedFrom = edge[1];
            contractedTo = edge[0];
        } else {
            hashMapSmaller = hashStart;
            hashMapBigger = hashEnd;
            contractedFrom = edge[0];
            contractedTo = edge[1];
        }

        Set<Integer> set = new HashSet<Integer>();

        // below loop is for updating adjacency list after each contraction and also adding the edge with height at least k to the force contracting set
        for (int key : hashMapSmaller.keySet()) {
            if (key != contractedTo && key != -1) {
                // hashmap hamsayehaye v
                HashMap<Integer, Integer> adjacencyMap = graphAdjacencyList.get(key);
                if (adjacencyMap.containsKey(contractedTo)) {
                    int weight = adjacencyMap.get(contractedFrom) + adjacencyMap.get(contractedTo);
                    adjacencyMap.put(contractedTo, weight);
                    // Add the edge to the set for force contraction if weight is at least k and the edge is not already in the set
                    // we should just check for the edges that the wight of them increase after contraction
                    if (weight >= k) {
                        int[] edgWithWeightAtLeastK = new int[2];
                        edgWithWeightAtLeastK[0] = key;
                        edgWithWeightAtLeastK[1] = contractedTo;
                        if (!isEdgeInSet(edgWithWeightAtLeastK)) {
                            contractionList.add(edgWithWeightAtLeastK);
                        }
                    }
                } else {
                    adjacencyMap.put(contractedTo, adjacencyMap.get(contractedFrom));
                }
                set.add(key);
            }
        }

        // Below code updates the edge lists that need to be force contracted if one of the end points are contracted from
        Iterator<int[]> itr = contractionList.iterator();
        while (itr.hasNext()) {
            int[] nextEdge = itr.next();
            if (nextEdge[0] == contractedFrom) {
                nextEdge[0] = contractedTo;
            }
            if (nextEdge[1] == contractedFrom) {
                nextEdge[1] = contractedTo;
            }
        }

        // remove edge in adjacency list when they are contracted
        int edgeWeight = hashStart.get(edge[1]);
        hashStart.remove(edge[1]);
        hashStart.put(-1, hashStart.get(-1) - edgeWeight);
        hashEnd.remove(edge[0]);
        hashEnd.put(-1, hashEnd.get(-1) - edgeWeight);

        // merge the smaller hash map into the bigger one
        for (HashMap.Entry<Integer, Integer> entry : hashMapSmaller.entrySet()) {
            int key = entry.getKey();

            if (key == -1) {
                continue;
            }

            //value is weight for smaller vertex
            int value = entry.getValue();
            if (hashMapBigger.containsKey(key)) {
                hashMapBigger.put(key, value + hashMapBigger.get(key));
            } else {
                hashMapBigger.put(key, value);
            }
            hashMapBigger.put(-1, hashMapBigger.get(-1) + value);
        }

        for (int index : set) {
            removeFromHashMap(graphAdjacencyList.get(index), contractedFrom);
        }

        graphAdjacencyList.remove(contractedFrom);

        if (calculateDegree(graphAdjacencyList, contractedTo) < k) {
            workSet.add(contractedTo);
        }
    }

    // Remove key from hash map
    public static void removeFromHashMap(HashMap map, int key) {
        Iterator<Map.Entry<Integer, Integer>> itr = map.entrySet().iterator();
        while (itr.hasNext()) {
            Map.Entry<Integer, Integer> entry = itr.next();
            if (entry.getKey() == key) {
                itr.remove();  // Call Iterator's remove method.
            }
        }
    }

    // This returns a random edge from given graph
    public static int[] pickRandomEdge(HashMap<Integer, HashMap<Integer, Integer>> graph, int h) {
        List<Integer> keysAsArray = new ArrayList(graph.keySet());
        List<Integer> innerKeysAsArray;
        Random r = new Random();
        int[] edge = new int[2];
        while (true) {
            edge[0] = keysAsArray.get(r.nextInt(keysAsArray.size()));

            HashMap<Integer, Integer> neighbors = graph.get(edge[0]);
            //dont choose -1
            if (neighbors.size() > 1) {
                innerKeysAsArray = new ArrayList<Integer>(neighbors.keySet());
                while (true) {
                    edge[1] = innerKeysAsArray.get(r.nextInt(innerKeysAsArray.size()));
                    if (edge[1] != -1) {
                        return edge;
                    }
                }
            }
        }
    }

    // This is the Algorithm 1 in paper
    public static Set<HashMap<Integer, HashMap<Integer, Integer>>> ContractAndCut(HashMap<Integer, HashMap<Integer, Integer>> G, int h) {
        //vaghti ke roye zirgraph sedash mianim bazam n verex darim dar sorati ke kamtar mishe!
        UF unionFind = new UF(n + 1);

        // This is a list of all the output components
        Set<HashMap<Integer, HashMap<Integer, Integer>>> allInducedSubGraphsOfRemovedVertices = new HashSet();

        // Copy the graph G into graph copyG
        HashMap<Integer, HashMap<Integer, Integer>> copyG = graphDeepCopy(G);
        while (!copyG.isEmpty()) {
            if (!workSet.isEmpty()) {

                while (!workSet.isEmpty()) {

                    int u = workSet.iterator().next();

                    //System.out.println("u" +u);
                    workSet.remove(u);

                    // Collect original vertices contracted to u

                    Set<Integer> originalVerticesContractedTou = unionFind.collectComponent(unionFind.find(u));
                    // if(originalVerticesContractedTou.size()>1){

                    // Build the induced sub graph of G[U]
                    HashMap<Integer, HashMap<Integer, Integer>> inducedSubGraph = buildInducedSubGraph(originalVerticesContractedTou);
                    inducedSubGraph = recursivelyRemoveVertexWithDegreeLessK(inducedSubGraph, k);

                    // Add the induced sub graph of G[U]  to the list of all graphs
                    allInducedSubGraphsOfRemovedVertices.add(inducedSubGraph);//}

                    //output(inducedSubGraph);

                    //let's delete u from all it's neighbors adjacency lists
                    Object[] neighbors = copyG.get(u).keySet().toArray();
                    for (int i = 0; i < neighbors.length; i++) {
                        int v = (int) neighbors[i]; //neighbor of u
                        if (v == -1) {
                            continue;
                        }
                        if (v == u) //self-loop
                        {
                            continue;
                        }

                        //delete u from adjacency list of v and reduce the degree of v
                        HashMap<Integer, Integer> hashMapOfV = copyG.get(v);
                        int weight = hashMapOfV.get(u);
                        hashMapOfV.remove(u);
                        hashMapOfV.put(-1, hashMapOfV.get(-1) - weight);

                        //if degree of v is smaller than k, add v to queue
                        if (calculateDegree(copyG, v) < k) {
                            workSet.add(v);
                        }
                    }

                    //now let's remove u from G
                    copyG.remove(u);
                }
            } else {

                if (copyG.size() == 1) {
                    copyG.clear();
                    continue;
                }
                int[] randomEdge = pickRandomEdge(copyG, h);

                //System.out.println("random edge:" + randomEdge[0] + randomEdge[1]);
                forcedContract(copyG, randomEdge);
                unionFind.union(randomEdge[0], randomEdge[1]);

                // Below code is doing the force contraction and contract any edges with weight at lease K
                while (!contractionList.isEmpty()) {

                    int[] e = contractionList.poll();
                    if (e[0] == e[1]) {
                        continue;
                    }
                    forcedContract(copyG, e);
                    unionFind.union(e[0], e[1]);
                }

                // System.out.println("afte force");
            }
        }
        return allInducedSubGraphsOfRemovedVertices;
    }

    public static void Decompose(HashMap<Integer, HashMap<Integer, Integer>> G) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.currentTimeMillis();
        ArrayList<Set<HashMap<Integer, HashMap<Integer, Integer>>>> pi = new ArrayList<>();
        Set piZero = new HashSet();
        piZero.add(G);
        pi.add(piZero);
        Set<HashMap<Integer, HashMap<Integer, Integer>>> piI;
        int components = 0;
        int new_components = 0;
        int h = 1;
        while (h <= t) {

            System.out.println(h + " itteration starts");

            pi.add(h, new HashSet());//majmoe khali to pi
            piI = pi.get(h);

            for (HashMap<Integer, HashMap<Integer, Integer>> Gprim : pi.get(h - 1)) {
                components++;
                piI.addAll(ContractAndCut(Gprim, h));
            }
            new_components = piI.size();

            h++;

        }
        long stopTime = System.currentTimeMillis();

        //writer = new PrintWriter("result.txt");
        PrintWriter writer2 = new PrintWriter("the-file-name2.txt", "UTF-8");

        Set<HashMap<Integer, HashMap<Integer, Integer>>> piIteration = pi.get(pi.size() - 1);
        int counter = 0;
        for (HashMap<Integer, HashMap<Integer, Integer>> graph : piIteration) {
            // output(graph);
            for (Map.Entry<Integer, HashMap<Integer, Integer>> entry : graph.entrySet()) {
                int vertex = entry.getKey();
                writer2.print(vertex);
                writer2.print(" ");

                System.out.print(vertex + "  ");
            }
            writer2.println();
            System.out.println();
            System.out.println();
            counter++;
            //System.out.println("Number of vertices in above sub graph: " + graph.size());
        }

        writer2.close();
        System.out.println("components:" + counter);
        System.out.print("Execution time: ");
        System.out.println((stopTime - startTime) / 1000.0);

    }

    public static void main(String[] argv) throws Exception {
        // Check if the correct number of arguments are provided
        if (argv.length < 3) {
            System.out.println("Usage: <program-name> <file-path> <k> <t>");
            System.exit(1);
        }

        // Parse the arguments
        String filePath = argv[0];
        k = Integer.parseInt(argv[1]);
        t = Integer.parseInt(argv[2]);

        // Create an ECC instance with the provided file path
        ECC ecc = new ECC(filePath);

        ecc.Decompose(mainAdjacencyList);
    }
}
