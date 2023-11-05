import java.util.Arrays;
import java.util.List;

public class ListLinearHeap {
    public int n; // number vertices
    public int key_cap; // the maximum allowed key value

    public int max_key; // possible max key
    public int min_key; // possible min key

    public int keys[]; // keys of vertices
    // keys[i] > key_cap if vertex i is not in the data structure

    public int heads[]; // head of doubly-linked list for a specific weight
    public int pres[]; // pre for doubly-linked list
    public int nexts[]; // next for doubly-linked list

    public ListLinearHeap(int _n, int _key_cap) {
        this.n = _n;
        this.key_cap = _key_cap;

        min_key = key_cap;
        max_key = 0;
        heads = keys = pres = nexts = null;
    }

    // initialize the data structure by (id, key) pairs
    // _n is the number of pairs, _key_cap is the maximum possible key value
    public void init(int _n, int _key_cap, int _ids[], int _keys[]) {
        if (keys == null) {
            keys = new int[n];
        }
        if (pres == null) {
            pres = new int[n];
        }
        if (nexts == null) {
            nexts = new int[n];
        }
        if (heads == null) {
            heads = new int[key_cap + 1];
        }
        min_key = _key_cap;
        max_key = 0;
        for (int i = 0; i <= _key_cap; i++) {
            heads[i] = n;
        }

        for (int i = 0; i < _n; i++) {
            insert(_ids[i], _keys[i]);
        }
    }

    // insert (id, key) pair into the data structure
    public void insert(int id, int key) {
        keys[id] = key;
        pres[id] = n;
        nexts[id] = heads[key];
        if (heads[key] != n) {
            pres[heads[key]] = id;
        }
        heads[key] = id;

        if (key < min_key) {
            min_key = key;
        }
        if (key > max_key) {
            max_key = key;
        }
    }

    // remove a vertex from the data structure
    public int remove(int id) {
        assert (keys[id] <= max_key);
        if (pres[id] == n) {
            assert (heads[keys[id]] == id);
            heads[keys[id]] = nexts[id];
            if (nexts[id] != n) {
                pres[nexts[id]] = n;
            }
        } else {
            int pid = pres[id];
            nexts[pid] = nexts[id];
            if (nexts[id] != n) {
                pres[nexts[id]] = pid;
            }
        }

        return keys[id];
    }

    int get_n() {
        return n;
    }

    int get_key_cap() {
        return key_cap;
    }

    int get_key(int id) {
        return keys[id];
    }

    public void get_ids(List<Integer> ids) {
        ids.clear();
        tighten();
        for (int i = min_key; i <= max_key; i++) {
            for (int id = heads[i]; id != n; id = nexts[id]) {
                ids.add(id);

            }
        }
    }

    public void get_ids_keys(List<Integer> ids, List<Integer> _keys) {
        ids.clear();
        _keys.clear();
        tighten();
        for (int i = min_key; i <= max_key; i++) {
            for (int id = heads[i]; id != n; id = nexts[id]) {
                ids.add(id);
                _keys.add(id);
            }
        }
    }

    public boolean empty() {
        tighten();
        return min_key > max_key;
    }

    public int size() {
        tighten();
        int res = 0;
        for (int i = min_key; i <= max_key; i++) {
            for (int id = heads[i]; id != n; id = nexts[id]) {
                ++res;
            }
        }
        return res;
    }

    // get the (id,key) pair with the maximum key value; return true if success, return false otherwise
    public boolean get_max(WrapInt id, WrapInt key) {
        if (empty()) {
            return false;
        }

        id.value = heads[max_key];
        key.value = max_key;
        return true;
    }

    // pop the (id,key) pair with the maximum key value; return true if success, return false otherwise
    public boolean pop_max(WrapInt id, WrapInt key) {
        if (empty()) {
            return false;
        }

        id.value = heads[max_key];
        key.value = max_key;
        heads[max_key] = nexts[id.value];
        if (heads[max_key] != n) {
            pres[heads[max_key]] = n;
        }
        return true;
    }

    // get the (id,key) pair with the minimum key value; return true if success, return false otherwise
    boolean get_min(WrapInt id, WrapInt key) {
        if (empty()) {
            return false;
        }

        id.value = heads[min_key];
        key.value = min_key;
        return true;
    }

    // pop the (id,key) pair with the minimum key value; return true if success, return false otherwise
    public boolean pop_min(WrapInt id, WrapInt key) {
        if (empty()) {
            return false;
        }

        id.value = heads[min_key];
        key.value = min_key;

        heads[min_key] = nexts[id.value];
        if (heads[min_key] != n) {
            pres[heads[min_key]] = n;
        }
        return true;
    }

    // increment the key of vertex id by inc
    public int increment(int id, int inc) {
        inc = 1;
        assert (keys[id] + inc <= key_cap);

        int new_key = keys[id] + inc;

        remove(id);
        insert(id, new_key);

        return new_key;
    }

    // decrement the key of vertex id by dec
    int decrement(int id, int dec) {
        dec = 1;
        assert (keys[id] >= dec);

        int new_key = keys[id] - dec;

        remove(id);
        insert(id, new_key);

        return new_key;
    }

    private void tighten() {
        while (min_key <= max_key && heads[min_key] == n) {
            ++min_key;
        }
        while (min_key <= max_key && heads[max_key] == n) {
            --max_key;
        }
    }
};
