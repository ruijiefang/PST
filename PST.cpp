/* Program Structure Tree */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cassert>
#include <string>
#include <set>
#include <cmath>
#include <vector>
#include <map>

struct Vertex {
  size_t ID;
  std::vector<size_t> Out;
};

struct Digraph {
  int N;
  Vertex* entryBlock;
  std::vector<Vertex*> Vertices;

  Vertex* getEntryBlock() { return this->entryBlock; }
  void    setEntryBlock(Vertex* entry) { this->entryBlock = entry; }
  void    Print() {
    using namespace std;
    cout << "+------------------- graph of " << N<< " nodes\n";
    for (size_t I=0; I<Vertices.size(); ++I) {
      cout << "Vertex " << Vertices[I]->ID << "\n";
      cout << " Neighbors: ";
      for(size_t J : Vertices[I]->Out) {
        Vertex *Q = Vertices[J];
        cout << Q->ID;
        cout << " ";
      }
      cout << "\n";
    }
    cout << "-------------------+\n";
  }

  ~Digraph() {
    for (size_t I = 0; I < Vertices.size(); ++I) {
      delete Vertices[I];
    }
  }

};

struct Bracket;

// Edge structure for the undirected
// multigraph version of control-flow graph
// required by the cycle-equivalence algorithm.
struct MultiEdge {
  // Source vertex
  unsigned U;

  // Sink vertex
  unsigned V;

  // Unique edge ID
  unsigned ID;

  // Equivalence class of edge
  unsigned Class;

  // most recent size of stack which
  // this edge is on top.
  unsigned RecentSize;

  // most recent equivalence class.
  unsigned RecentClass;

  // Is this edge a reversed edge
  bool Rev;

  // Is this edge a capping edge
  bool IsCappingEdge;

  // A reverse pointer into the
  // container carrying this edge inside
  // the linked list of back-edges.
  // There is a one-to-one correspondence
  // between a MultiEdge and a Bracket.
  Bracket *B;

  void Init() {
    this->Class = 0;
    this->RecentSize = 0;
    this->RecentClass = 0;
    this->IsCappingEdge = false;
    this->B = nullptr;
  }

  MultiEdge() {
    U = 0; V = 0; Rev = false; ID = 0;
    Init();
  }

  MultiEdge(unsigned U, unsigned V) {
    this->U = U;
    this->V = V;
    this->Rev = false;
    this->ID = 0;
    Init();
  }

  MultiEdge(unsigned U, unsigned V, unsigned ID) {
    this->U = U;
    this->V = V;
    this->Rev = false;
    this->ID = ID;
    Init();
  }

  bool operator==(const MultiEdge& RHS) {
    return this->ID == RHS.ID;
  }

  Bracket* getBracket() const {
    assert(this->B != nullptr && "Getting a null bracket from a MultiEdge.");
    return this->B;
  }

  void setBracket(Bracket* B) {
    this->B = B;
  }
};

// A container for backedges inside the
// linked list (BracketList) of back-edges maintained
// by the cycle equivalence algorithm.
struct Bracket {

  // A pointer to the edge we're carrying.
  MultiEdge *E;

  // Next in list.
  Bracket *Next;

  // Previous in list.
  Bracket *Prev;

  Bracket() {
    Next = Prev = nullptr;
  }

  Bracket (MultiEdge *E) {
    this->E = E;
    Next = Prev = nullptr;
  }

  Bracket (MultiEdge *E, Bracket *Next, Bracket *Prev) {
    this->E = E;
    this->Next = Next;
    this->Prev = Prev;
  }
};

// A linked list structure of back-edges, required
// by the cycle equivalence algorithm supporting O(1) deletion.
// All back-edges from a vertex u's descendants that reach over u
// are called brackets.
struct BracketList {

  // The topmost element in list.
  Bracket * Top;

  // The bottom-most element in list.
  Bracket * Bottom;

  // The size of the current list.
  size_t Size;

  BracketList() {
    this->Top = nullptr;
    this->Bottom = nullptr;
    this->Size = 0;
  }

  // Polymorphic Push() method. Pushes
  // a heap-allocated bracket onto the list.
  void Push(Bracket *B) {
    Bracket *T = this->Top;
    this->Top = B;
    B->Prev = nullptr;
    B->Next = T;
    T->Prev = B;
    if (Bottom == nullptr)
      Bottom = B;
  }

  // Polymorphically receive a MultiEdge
  // reference; build its container Bracket
  // on the heap and insert it into the list.
  // A deleted bracket is always de-allocated
  // automatically in the Delete() method.
  void Push(MultiEdge& E) {
    Bracket *B = new Bracket(&E);
    this->Push(B);
    E.setBracket(B);
  }

  // Deletes a bracket pointer B in O(1)-time
  // by unlinking it from the list.
  // Also deallocates its memory since all
  // bracket structures are heap-allocated.
  void Delete(Bracket *B) {
    if (B->Prev == nullptr) {
      // B is Top node.
      this->Top = B->Next;
      B->Next->Prev = nullptr;
      B->Prev = B->Next = nullptr;
    } else if (B->Next == nullptr) {
      // B is Bottom node.
      this->Bottom = B->Prev;
      B->Prev->Next = nullptr;
      B->Prev = B->Next = nullptr;
    } else {
      // B is somewhere in the middle.
      B->Next->Prev = B->Prev;
      B->Prev->Next = B->Next;
      B->Prev = B->Next = nullptr;
    }
    delete B;
  }

  // Concatenates a sublist into this current list.
  // After concatenation, the sublist structures are
  // permanently damaged and cannot be reused.
  void Concat(const BracketList& NextList) {
    // Trivial case.
    if (NextList.Size == 0) return;
    this->Size += NextList.Size;
    this->Bottom->Next = NextList.Top;
    NextList.Top->Prev = this->Bottom;
    this->Bottom = NextList.Bottom;
  }
};

// Use DFS to calculate DFS numbers, parent pointers, and back-edges for
// the cycle-equivalence algorithm.
void UndirectedDepthFirstTraversal(const std::vector< std::vector<MultiEdge> >& UGraph,
                                   std::vector< std::vector<MultiEdge> >& BackEdgesInto,
                                   std::vector<unsigned>& Parent,
                                   std::vector<unsigned>& DFSNum,
                                   std::vector<unsigned>& DFSOrder,
                                   unsigned& DFSCounter, int U) {
  // Visited?
  if (DFSNum[U] != 0)
    return;

  DFSNum[U] = DFSCounter++;
  DFSOrder.push_back(U);

  // Visit our neighbors.
  for (const MultiEdge &E : UGraph[U]) {
    if (DFSNum[E.V] != 0) {
      // Back edge found. Mark it.
      BackEdgesInto[E.V].push_back(E);
    }
    // Otherwise, visit the unvisited edge.
    Parent[E.V] = E.U;
    UndirectedDepthFirstTraversal(UGraph, BackEdgesInto, Parent, DFSNum, DFSOrder, DFSCounter, E.V);
  }
}

// O(|E|)-time canonical SESE region calculation using cycle-equivalence algorithm.
std::vector< std::set<Vertex*> > calculateCanonicalSESERegions(Digraph &G) {
  // Establish a bijective mapping between BBs in the CFG to natural numbers.
  std::map<Vertex *, unsigned> VertexNumbering;
  std::vector<Vertex *> VertexNumberToBB;
  
  // Use the numbering to calculate directed/undirected adjacency lists
  // for the cycle equivalence algorithm.
  std::vector< std::vector<MultiEdge> > DirectedAdjList;
  std::vector< std::vector<MultiEdge> > UndirectedAdjList;
  std::vector<MultiEdge> DirectedEdgeList;
  unsigned VertexNumber = 1;
  unsigned EdgeNumber = 1;

  // First, map each BB to numbers.
  for (Vertex* BB : G.Vertices)
    VertexNumbering[BB] = VertexNumber++;

  // Initialize adjacency list.
  DirectedAdjList.assign(VertexNumber + 1, std::vector<MultiEdge>());
  UndirectedAdjList.assign(VertexNumber + 1, std::vector<MultiEdge>());

  // Next, populate the adjacency list using CFG information.
  for(Vertex *BB : G.Vertices) {
    unsigned SrcVertex = VertexNumbering[BB];
    VertexNumberToBB.push_back(BB);
    for (unsigned I = 0; I < BB->Out.size(); ++I) {
      Vertex *Succ = G.Vertices[BB->Out[I]];
      unsigned DstVertex = VertexNumbering[Succ];
      DirectedAdjList[SrcVertex].push_back(MultiEdge(SrcVertex, DstVertex, EdgeNumber));
      DirectedEdgeList.push_back(MultiEdge(SrcVertex, DstVertex, EdgeNumber++));
    }
  }

  // Establish source, sink nodes of the CFG.
  // If the CFG itself is SEME, use 0 as a "virtual" terminating block.
  unsigned Source = VertexNumbering[G.getEntryBlock()], Sink;
  unsigned NumExits = 0;
  for (unsigned I = 1; I < DirectedAdjList.size(); ++I)
    if (DirectedAdjList[I].size() == 0){
      ++NumExits;
      Sink = I;
    }

  assert(NumExits > 0 && "Number of control-flow graph exits must > 0.");

  if (NumExits > 1) {
    Sink = 0;
    // Link all exit nodes to 0.
    for (unsigned I = 1; I < DirectedAdjList.size(); ++I)
      if (DirectedAdjList[I].size() == 0){
        DirectedAdjList[I].push_back(MultiEdge(I, Sink, EdgeNumber));
        DirectedEdgeList.push_back(MultiEdge(I, Sink, EdgeNumber++));
      }
  }

  // Now convert the directed graph into an undirected (multi)-graph for finding
  // SESE regions, and add edge from Sink to Source.
  for (MultiEdge& E : DirectedEdgeList) {
    MultiEdge ReverseE = MultiEdge(E.V, E.U, E.ID);
    ReverseE.Rev = true;
    UndirectedAdjList[E.U].push_back(E);
    UndirectedAdjList[E.V].push_back(ReverseE);
  }
  MultiEdge STEdge(Source, Sink, 0);
  MultiEdge TSEdge(Sink, Source, 0);
  TSEdge.Rev = true;
  UndirectedAdjList[Source].push_back(STEdge);
  UndirectedAdjList[Sink].push_back(TSEdge);

  // Set up vertex attributes maintained by our algorithm.
  std::vector<unsigned> UndirectedDFSNum(UndirectedAdjList.size());
  std::vector<BracketList> BList(UndirectedAdjList.size());
  std::vector<unsigned> Hi(UndirectedAdjList.size());
  std::vector<unsigned> Parent(UndirectedAdjList.size());

  // Set up information maintained by DFS traversal of graph.
  std::vector<unsigned> DFSNum(UndirectedAdjList.size());
  std::vector<unsigned> DFSOrder;
  std::vector< std::vector<MultiEdge> > BackEdgesInto(UndirectedAdjList.size());
  // Counter for assigning DFS number to each vertex.
  unsigned DFSCounter = 1;
  // Counter for registering equivalence classes of SESE regions.
  unsigned EquivalenceClassCounter = 1;

  // Do a DFS traversal of the graph to set up
  Parent[Source] = 0;
  UndirectedDepthFirstTraversal(UndirectedAdjList, BackEdgesInto, Parent, DFSNum, DFSOrder, DFSCounter, Source);

  // Calculate SESE regions in reverse depth-first order.
  std::reverse(DFSOrder.begin(), DFSOrder.end());
  for (const unsigned& N : DFSOrder) {
    // Calculate Hi0-2 and HiChild for capping back-edges.
    unsigned Hi0 = 0xFFFFFFFEU, Hi1 = 0xFFFFFFFEU, HiChild, Hi2 = 0xFFFFFFFEU;
    for (const MultiEdge& E : UndirectedAdjList[N]) {
      if (DFSNum[E.V] > DFSNum[E.U] && Hi0 >= DFSNum[E.V]) {
        // Update Hi0 := min{dfsnum(t); (n, t) backedge}
        Hi0 = DFSNum[E.V];
      } else if (DFSNum[E.V] < DFSNum[E.U] && Hi1 >= Hi[E.V]) {
        // V is child of U.
        // This case is always triggered inductively. The base case
        // will never rest in this if-condition, since the base
        // case is a leaf in the DFS spanning tree.
        // Here, update Hi1 = min{c.Hi; c child of N}
        Hi1 = Hi[E.V];
      }
    }

    // Hi of N = min{Hi0, Hi1}
    if (Hi0 < Hi1)
      Hi[N] = Hi0;
    else
      Hi[N] = Hi1;

    // Compute HiChild = any child c of N having Hi[c] = Hi[N]
    for(const MultiEdge& E : UndirectedAdjList[N])
      if (DFSNum[E.V] > DFSNum[E.U] && Hi[E.V] == Hi[E.U]) {
        HiChild = E.V;
        break;
      }

    // Compute Hi2 = min{C.Hi; C is child of N other than HiChild}
    for(const MultiEdge& E : UndirectedAdjList[N])
      if (DFSNum[E.V] > DFSNum[E.U] && E.V != HiChild && Hi2 > Hi[E.V])
        Hi2 = Hi[E.V];

    // Next, compute the BracketList of node N.
    // First, concatenate together all bracket lists of children of N
    // into the bracket list of the current node.
    // Note that, the bracket lists of children will be destroyed
    // and no longer valid after we consolidated everything into the
    // current bracket list.
    for(const MultiEdge& E : UndirectedAdjList[N])
      if (DFSNum[E.V] > DFSNum[E.U])
        BList[N].Concat(BList[E.V]);
    // Next, delete back-edges from descendants of N to N.
    // If a back-edge B is a capping back-edge, we simply delete
    // it. Otherwise, we assign a new equivalence class to it
    // after deleting it.
    for(MultiEdge& E : BackEdgesInto[N]) {
      BList[N].Delete(E.getBracket());
      E.setBracket(nullptr);
      if (!E.IsCappingEdge && E.Class == 0)
        E.Class = EquivalenceClassCounter++;
    }
    // For each backedge B out of N,
    // Push B into the BracketList of N.
    for(MultiEdge& E : UndirectedAdjList[N])
      if (DFSNum[E.V] > DFSNum[E.U]) {
        BList[N].Push(E);
      }
    // If Hi2 < Hi0 for the current node N,
    // then we must create a capping back-edge from
    // N to node Hi2.
    if (Hi2 < Hi0) {
      MultiEdge CappingEdge(N, Hi2);
      CappingEdge.IsCappingEdge = true;
      BList[N].Push(CappingEdge);
    }

    // Finally, determine the equivalence class for the in-edge
    // from Parent[N] to node N (if N is not Source).
    if (Parent[N] != 0) {
      for (MultiEdge& E : UndirectedAdjList[N])
        if (E.V == Parent[N]) {
          Bracket * B = BList[N].Top;
          assert(B != nullptr && "Empty Bracket List");
          if (B->E->RecentSize != BList[N].Size) {
            B->E->RecentSize = BList[N].Size;
            B->E->RecentClass = EquivalenceClassCounter++;
          }
          E.Class = B->E->RecentClass;
          // Check for E, B equivalence:
          if (B->E->RecentSize == 1)
            B->E->Class = E.Class;
          break;
        }
    }
  }

  // After finished calculating cycle equivalence,
  // map equivalence classes of SESE regions back to pairs
  // of LLVM BasicBlock objects.

  std::vector< std::set<Vertex*> > EquivalentCycles(EquivalenceClassCounter);
  for (unsigned I = 1; I < UndirectedAdjList.size(); ++I) {
    for (const MultiEdge& E : UndirectedAdjList[I]) {
      unsigned EquivClass = E.Class;
      Vertex *SrcBlock = VertexNumberToBB[E.U], *SinkBlock = VertexNumberToBB[E.V];
      EquivalentCycles[EquivClass].insert(SrcBlock);
      EquivalentCycles[EquivClass].insert(SinkBlock);
    }
  }
  return EquivalentCycles;
}

// main entry.
int main(int argc, char **argv)
{
  using namespace std;
  cout << "PST.cpp "<<" taking input from stdin.\n";
  
  Digraph G;
  cin >> G.N;
  G.Vertices.assign(G.N, nullptr);
  for(size_t I = 0; I < G.N; ++I) {
    Vertex *V = new Vertex;
    cin >> V->ID;
    size_t NeighborCnt;
    cin >> NeighborCnt;
    for(size_t J = 0; J < NeighborCnt; ++J) {
      size_t Neighbor;
      cin >> Neighbor;
      V->Out.push_back(Neighbor);
    }
    G.Vertices[V->ID] = V;
    cout<< "Created vertex " << V->ID << "\n";
  }
  G.Print();
  cout<<"starting calculation...\n";
  vector< set<Vertex*> > Classes = calculateCanonicalSESERegions(G);
  cout<<"finished calculation.\n";
  cout<<"number of equivalences:"<<Classes.size()<<"\n";
  return 0;
}


