import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

public class UF {
    private int[] id;    // id[i] = parent of i
    private int[] sz;    // sz[i] = number of objects in subtree rooted at i
    private int count;   // number of components

    /**
     * Create an empty union find data structure with N isolated sets.
     */
    public UF(int N) {
        count = N;
        id = new int[N];
        sz = new int[N];
        for (int i = 0; i < N; i++) {
            id[i] = i;
            sz[i] = 1;
        }
    }

    /**
     * Return the id of component corresponding to object p.
     */
    public int find(int p) {
        while (p != id[p]) {
            id[p] = id[id[p]];
            p = id[p];
        }
        return p;
    }

    /**
     * Return the number of disjoint sets.
     */
    public int count() {
        return count;
    }

    /**
     * Are objects p and q in the same set?
     */
    public boolean connected(int p, int q) {
        return find(p) == find(q);
    }

    /**
     * Replace sets containing p and q with their union.
     */
    public void union(int p, int q) {
        int i = find(p);
        int j = find(q);
        if (i == j) {
            return;
        }

        // make smaller root point to larger one
        if (sz[i] < sz[j]) {
            id[i] = j;
            sz[j] += sz[i];
        } else {
            id[j] = i;
            sz[i] += sz[j];
        }
        count--;
    }

    //it returuns hashmap that keys are roots value are set of child of this roots and root itself
    public Map<Integer, Set<Integer>> components() {
        Map<Integer, Set<Integer>> res = new HashMap<>();
        for (int p = 0; p < id.length; p++) {
            int rootp = find(p);
            if (!res.containsKey(rootp)) {
                res.put(rootp, new HashSet<>());
            }
            res.get(rootp).add(p);
        }
        return res;
    }

    public List<Integer> collectComponent(int compid) {
        Set<Integer> res = new HashSet<>();
        List<Integer> res2 = new ArrayList<>();
        for (int p = 0; p < id.length; p++) {
            if (find(p) == compid) {
                if (!res.contains(p)) {
                    res.add(p);
                    res2.add(p);
                }
            }
        }
        return res2;
    }

    public static void main(String[] args) throws Exception {
        long startTime = System.currentTimeMillis();
        int N = 4_847_571;
        UF uf = new UF(N);

        // read in a sequence of pairs of integers (each in the range 0 to N-1),
        // calling find() for each pair: If the members of the pair are not already
        // call union() and print the pair.
        BufferedReader br = new BufferedReader(
                new FileReader("C:\\Users\\hanieh\\Downloads\\UF - Copy-20211122T072339Z-001\\UF - Copy\\soc-Epinions1.txt"));
        String line;
        while ((line = br.readLine()) != null) {
            String[] lineArr = line.split("\\s+");
            int p = Integer.parseInt(lineArr[0]);
            int q = Integer.parseInt(lineArr[1]);

            if (uf.connected(p, q)) {
                continue;
            }

            uf.union(p, q);
        }

        System.out.println("# components: " + uf.count());
        br.close();

        long stopTime = System.currentTimeMillis();
        System.out.println((stopTime - startTime) / 1000.0);

        List<Integer> comp = uf.collectComponent(0);
        System.out.println(comp.size());

        stopTime = System.currentTimeMillis();
        System.out.println((stopTime - startTime) / 1000.0);
    }

}
