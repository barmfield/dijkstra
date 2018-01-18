
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.Collections;
import java.util.ListIterator;

public class Dijkstra {

    // Keep a fast index to nodes in the map
    private Map<String, Vertex> vertexNames;

    /**
     * Construct an empty Dijkstra with a map. The map's key is the name of a
     * vertex and the map's value is the vertex object.
     */
    public Dijkstra() {
        vertexNames = new HashMap<String, Vertex>();

    }

    /**
     * Adds a vertex to the dijkstra. Throws IllegalArgumentException if two
     * vertices with the same name are added.
     *
     * @param v (Vertex) vertex to be added to the dijkstra
     */
    public void addVertex(Vertex v) {
        if (vertexNames.containsKey(v.name)) {
            throw new IllegalArgumentException("Cannot create new vertex with existing name.");
        }
        vertexNames.put(v.name, v);
    }

    /**
     * Gets a collection of all the vertices in the dijkstra
     *
     * @return (Collection<Vertex>) collection of all the vertices in the
     * dijkstra
     */
    public Collection<Vertex> getVertices() {
        return vertexNames.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name (String) name of the vertex object requested
     * @return (Vertex) vertex object associated with the name
     */
    public Vertex getVertex(String name) {
        return vertexNames.get(name);
    }

    /**
     * Adds a directed edge from vertex u to vertex v
     *
     * @param nameU (String) name of vertex u
     * @param nameV (String) name of vertex v
     * @param cost (double) cost of the edge between vertex u and v
     */
    public void addEdge(String nameU, String nameV, Double cost) {
        if (!vertexNames.containsKey(nameU)) {
            throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
        }
        if (!vertexNames.containsKey(nameV)) {
            throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
        }
        Vertex sourceVertex = vertexNames.get(nameU);
        Vertex targetVertex = vertexNames.get(nameV);
        Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
        sourceVertex.addEdge(newEdge);
    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a
     * directed edge from u to v, then a directed edge from v to u
     *
     * @param nameU (String) name of vertex u
     * @param nameV (String) name of vertex v
     * @param cost (double) cost of the edge between vertex u and v
     */
    public void addUndirectedEdge(String nameU, String nameV, double cost) {
        addEdge(nameU, nameV, cost);
        addEdge(nameV, nameU, cost);
    }

    // STUDENT CODE STARTS HERE
    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux (double) x coordinate of point u
     * @param uy (double) y coordinate of point u
     * @param vx (double) x coordinate of point v
     * @param vy (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {

        return sqrt(pow((vx - ux), 2) + pow((vy - uy), 2));
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     *
     * line reading code copied from below.
     */
    public void computeAllEuclideanDistances() throws FileNotFoundException, IOException {
        String edgeFile = "citypairs.txt";
        BufferedReader edgeFileBr = new BufferedReader(new FileReader(edgeFile));
        String line;
        while ((line = edgeFileBr.readLine()) != null) {
            String[] parts = line.split(",");
            if (parts.length != 3) {
                edgeFileBr.close();
                throw new IOException("Invalid line in edge file " + line);
            }

            // get values for euclidean distance computation
            Vertex v1 = vertexNames.get(parts[0]);
            Vertex v2 = vertexNames.get(parts[1]);

            // Iterate through all edges in the linked list attached to each
            // vertex. If the source and target from each edge match the vertices
            // from the line reader above, update the distance on the edge
            // going from A to B.
            for (int i = 0; i < v1.adjacentEdges.size(); i++) {
                if (v1.adjacentEdges.get(i).source == v1
                        && v1.adjacentEdges.get(i).target == v2) {
                    v1.adjacentEdges.get(i).distance
                            = computeEuclideanDistance(v1.x, v1.y, v2.x, v2.y);
                }
            }

            // update distance on edge going from B to A
            for (int i = 0; i < v2.adjacentEdges.size(); i++) {
                if (v2.adjacentEdges.get(i).source == v2
                        && v2.adjacentEdges.get(i).target == v1) {
                    v2.adjacentEdges.get(i).distance
                            = computeEuclideanDistance(v1.x, v1.y, v2.x, v2.y);
                }
            }
        }
        edgeFileBr.close();
    }

    /**
     * Dijkstra's Algorithm.
     *
     * @param s (String) starting city name
     */
    public void doDijkstra(String s) {
        //add all unvisited vertices to a list
        LinkedList<Vertex> unvisited = new LinkedList(vertexNames.values());
        LinkedList<Vertex> visited = new LinkedList();
        // For all vertices, set distance to infiniti and known to false
        int size = unvisited.size();
        for (int i = 0; i < size; i++) {
            unvisited.get(i).distance = Double.POSITIVE_INFINITY;
            unvisited.get(i).known = false;
        }
        
        vertexNames.get(s).distance = 0;
        vertexNames.get(s).prev = null;
        while(!unvisited.isEmpty()) {            
            // find vertex with minimum distance from origin point
            double min = Integer.MAX_VALUE;
            Vertex v = null;
            for (int i = 0; i < unvisited.size(); i++) {
                if(unvisited.get(i).distance < min) {
                    v = unvisited.get(i); 
                    min = v.distance;
                }
            }
            
            v.known = true;

            // for each vertex adjacint to v,
                // if unknown: set distance from v
            for (int i = 0; i < v.adjacentEdges.size(); i++) {
                // variable to house adjacent vertices being cycled through
                Vertex w = v.adjacentEdges.get(i).target;
                
                if(!w.known){
                    double costVW = v.adjacentEdges.get(i).distance;  
                    // if v.distance plus the cost of the edge from v to w
                    // is less than w.distance, w.distance = v.distance + costVW
                    // update w's path
                    if(v.distance + costVW < w.distance){
                        w.distance = v.distance + costVW;
                        w.prev = v;
                    }
                }
            }
            // add to and remove from visited and unvisited lists
            visited.add(v);
            unvisited.remove(v);   
        }                   
    }

    /**
     * Returns a list of edges for a path from city s to city t. This will be
     * the shortest path from s to t as prescribed by Dijkstra's algorithm
     *
     * @param s (String) starting city name
     * @param t (String) ending city name
     * @return (List<Edge>) list of edges from s to t
     */
    public List<Edge> getDijkstraPath(String s, String t) {
        doDijkstra(s);        
        LinkedList<Edge> pathList = new LinkedList();
        
        
        Vertex last = vertexNames.get(t);
        
        // start with the destination vertex
        // using destination.prev, search adjacentEdges for edge connecting the two
        // add to list
        // replace current vertex with destination.prev
        // run until .prev == null
        
        while(last.prev != null) {
            Edge newEdge = new Edge(last.prev, last, 
                    computeEuclideanDistance(last.prev.x, last.prev.y, last.x, last.y));
            pathList.add(newEdge);
            last = last.prev;
        }
            
        Collections.reverse(pathList);
        return pathList;
    }

    // STUDENT CODE ENDS HERE
    /**
     * Prints out the adjacency list of the dijkstra for debugging
     */
    public void printAdjacencyList() {
        for (String u : vertexNames.keySet()) {
            StringBuilder sb = new StringBuilder();
            sb.append(u);
            sb.append(" -> [ ");
            for (Edge e : vertexNames.get(u).adjacentEdges) {
                sb.append(e.target.name);
                sb.append("(");
                sb.append(e.distance);
                sb.append(") ");
            }
            sb.append("]");
            System.out.println(sb.toString());
        }
    }

    /**
     * A main method that illustrates how the GUI uses Dijkstra.java to read a
     * map and represent it as a graph. You can modify this method to test your
     * code on the command line.
     */
    public static void main(String[] argv) throws IOException {
        String vertexFile = "cityxy.txt";
        String edgeFile = "citypairs.txt";

        Dijkstra dijkstra = new Dijkstra();
        String line;

        // Read in the vertices
        BufferedReader vertexFileBr = new BufferedReader(new FileReader(vertexFile));
        while ((line = vertexFileBr.readLine()) != null) {
            String[] parts = line.split(",");
            if (parts.length != 3) {
                vertexFileBr.close();
                throw new IOException("Invalid line in vertex file " + line);
            }
            String cityname = parts[0];
            int x = Integer.valueOf(parts[1]);
            int y = Integer.valueOf(parts[2]);
            Vertex vertex = new Vertex(cityname, x, y);
            dijkstra.addVertex(vertex);
        }
        vertexFileBr.close();

        BufferedReader edgeFileBr = new BufferedReader(new FileReader(edgeFile));
        while ((line = edgeFileBr.readLine()) != null) {
            String[] parts = line.split(",");
            if (parts.length != 3) {
                edgeFileBr.close();
                throw new IOException("Invalid line in edge file " + line);
            }
            dijkstra.addUndirectedEdge(parts[0], parts[1], Double.parseDouble(parts[2]));
        }
        edgeFileBr.close();

        // Compute distances. 
        // This is what happens when you click on the "Compute All Euclidean Distances" button.
        dijkstra.computeAllEuclideanDistances();

        // print out an adjacency list representation of the graph
        dijkstra.printAdjacencyList();

        // This is what happens when you click on the "Draw Dijkstra's Path" button.
        // In the GUI, these are set through the drop-down menus.
        String startCity = "Denver";
        String endCity = "Miami";

        // Get weighted shortest path between start and end city. 
        System.out.println("Getting Dijkstra's path...");
        List<Edge> path = dijkstra.getDijkstraPath(startCity, endCity);
        System.out.println("Path Completed.");
        
        System.out.print("Shortest path between " + startCity + " and " + endCity + ": ");
        System.out.println(path);
        
    }

}
