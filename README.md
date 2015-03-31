# mincost.cc

## Synopsis

mincost.cc is a C++ program which solves the minimum-cost flow problem on a discrete directed network given integer arc capacities, unit costs, 
and divergence values using the epsilon relaxation algorithm with epsilon scaling. Please see [wikipedia](http://en.wikipedia.org/wiki/Minimum-cost_flow_problem) for details on the problem and background. In simple terms, it finds a flow vector x indexed by the arcs of the graph that yields the minimum cost and 
which has the specified divergence values at every node.

### Notation:

1. `x` : flow vector indexed by the arcs of the graph. Each element of `x` corresponds to the flow assigned to a given arc
2. `d` : cost vector indexed by the arcs of the graph. Each element of `d` gives the unit cost for a given arc.
3. `c` : capacity vector indexed by arc. Each element of `c` gives the upper capacity of each arc. The lower capacities are all assumed to be 0.
4. `b` : supply vector indexed by node. Each element of `b` gives the **desired supply** at each node of the graph.
5. `r` : reduced cost vector indexed by arc. `r[j] = d[j] + u[i] - u[i']` where `j ~ (i, i')` is the arc from node `i` to node `i'`.
6. `s` : surplus vector indexed by node. Indicates the surplus supply at each node. i.e. 
  ```s[i] = (net flow into node  i) - b[i].``` 
7. `u` : node potential vector indexed by node.

### Data Format

1. Row 1 indicates |N| and |A|, where N is the node set and A is the arc set.
2. Rows 2 to |A|+1 correspond to arcs, indicating the start node, end node, cost, and upper capacity of each arc (lower capacity is assumed to be zero). 
3. The last |N| rows correspond to nodes, indicating the supply at each node.

There are four data sets supplied with this repository, with the following properties:

1. mcnf1.dat: |N|=4, |A|=5, optimal cost = 14
2. mcnf2.dat: |N|=20, |A|=80, optimal cost = 30464
3. mcnf3.dat: |N|=49, |A|=520, optimal cost = 1.73566449E+11 (corrected) (caution: this network has some parallel arcs)
4. mcnf4.dat: |N|=10000, |A|=60000, optimal cost = 21514702 (maybe)

## Algorithm Description:

### Complementary Epsilon Slackness Condition:

A given flow `x` and potential `u` satisfy epsilon-complementary slackness if for all arcs j ~ (i,i'), `d[j] + u[i] - u[i']` is

* in the interval `[-\epsilon, \epsilon]` if `0 < x[j] < c[j]`
* `>= -\epsilon` if `0 = x[j]`
* `<= \epsilon` if `x[j] = c[j]`

### Initialization: 
1. Set `u[i] = 0` for all nodes `i`.
2. For each arc `j ~ (i,i')`, compute `r[j] = d[j] + u[i] - u[i']`
3. Initialize flow `x` by calling `initflow(x,r,c,A)`
4. Set epsilon = maximum degree of any node so that x and u satisfy epsilon-complementary slackness condition.
5. Compute node surplus `s[i] = (net flow into node i) - b[i]`

### Main algorithm:
1. If `s[i] = 0` for all nodes `i`, then `x` is feasible, so stop. Else pick any node ibar with `s[ibar] > 0`. 
2. While `s[ibar] != 0`:
  1. Paint each arc `j ~ (i,i')`:
    * white if `-epsilon <= r[j] <= -epsilon/2` and `x[j] < c[j]`
    * black if `-epsilon/2 <= r[j] <= epsilon` and `x[j] > 0`
    * red otherwise
  2. If there is a black arc `jbar` out of `ibar`, then decrease `x[jbar]` by `min(x[jbar], s[ibar])`. Otherwise, if there is a white arc `jbar` into `ibar`, then increase `x[jbar]` by `min( c[jbar] - x[jbar], s[ibar])`.
  3. If no criteria from previous step apply, then increase `u[ibar]` by as much as possible while maintaining `epsilon`-complementary slackness. In particular, set `u[ibar] = u[ibar] + alpha`, where
  ```cplusplus
  alpha1 = min(-r[j] + epsilon) such that x[j] > 0 and j ~ (ibar, i')
  alpha2 = min(r[j] + epsilon) such that x[j] < c[j] and j ~ (i', ibar)
  alpha = min(alpha1, alpha2)
  ```
3. `epsilon = epsilon/2`
4. Repeat 1 until `epsilon < 1/N`.

## Code Example

```terminal
$ g++ mincost.cc -o mincost
$ ./mincost mcnf1.dat
Run time: 0.000361 seconds
Current flow is feasible
Minimum cost is: 14
$ ./mincost mcnf4.dat
Run time: 4.459991 seconds
Current flow is feasible
Minimum cost is: 21514702
```

## Author

- Mark Hubenthal
- markhubenthal@gmail.com
- [hubenjm.github.io](hubenjm.github.io)

