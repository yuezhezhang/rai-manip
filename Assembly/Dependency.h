#include <Core/array.h>

#include <KOMO/komo.h>
#include <Kin/proxy.h>
#include <Geo/fclInterface.h>

/*struct NeighborHandler{
  // TODO: remove
  //intA groundPartsIDs = {312, 2, 0, 131, 130, 129, 128, 127, 126, 125, 124, 123, 375, // Group 1
  //                       332, 40, 39, 176, 41, 43, 333, // Group 2
  //                       351, 77, 76, 214, 215, 216, 217, 78, 80, 352}; // Group 3
  intA groundPartsIDs = {1, 2, 3, 4, 18, 20, 17, 19};
  //intA groundPartsIDs = {1};

  StringA groundParts;
  intA neighbors;

  NeighborHandler(){
    for (auto g: groundPartsIDs){
      //groundParts.append(STRING("o"<<g));
      groundParts.append(STRING("b"<<g));
    }
  };

  StringA getNeighbours(const rai::String &objName) const{
    const uint id = std::stoi(objName.getSubString(1, -1).p);
    //const intA neighbor_ids = neighbors[id];
    const intA neighbor_ids = neighbors[id-1];

    StringA tmp;
    for (auto a: neighbor_ids){
      if (a >= 0){
        //tmp.append(STRING("o"<<a));
        tmp.append(STRING("b"<<a));
      }
    }

    return tmp;
  }

  rai::Array<StringA> getConnectedComponents(const StringA &feasible) const{
    rai::Array<StringA> cc;

    StringA fCpy = feasible;

    while (fCpy.N > 0){
      StringA toAdd = {fCpy.popFirst()};
      StringA component;

      while (toAdd.N > 0){
        const rai::String obj = toAdd.popFirst();
        component.append(obj);
        
        // add all neighbors to the same connected component
        const StringA n = getNeighbours(obj);
        for (const rai::String pn: n){
          if (!fCpy.contains(pn)) {continue;} // the part was already added
          if (component.contains(pn) || toAdd.contains(pn)) {continue;} // the part is already in the connected component

          toAdd.append(pn);
          fCpy.removeValue(pn);
        }
      }

      if (component.N == 1){std::cout << "A" << component << getNeighbours(component(0)) << std::endl;}
      cc.append(component);
    }

    return cc;
  }
};*/

struct NeighborHandler{
  rai::Configuration C;

  uint N_parts;                 // number of parts (including floor)
  intA neighborMatrix;          // adjacency Matrix (N_parts x N_parts)
  intA neighborList;            // list of neighbors for each part
  intA neccessaryNeighbors;     // A_ij == 1 if i and j are neighbors and i_z < j_z; A_ij == 2 if i and j are neighbors and j_z < i_z
  StringA shapeNames;           // list of part names in order of model file
  StringA groundParts;          // list of part names that touch the ground
  doubleA distanceMatrix;       // pairwise part distances (N_parts x N_parts)

  NeighborHandler(const rai::Configuration& _C){
    C.copy(_C);
    addFloor();

    for(rai::Frame* f:C.frames) if(f->shape && f->name!="floor" && !f->name.startsWith("_")){
      f->setJoint(rai::JT_rigid);
      f->setContact(1.);
      //f->setMass(.1);
    }

    C.stepSwift();
    setShapeNames();
    N_parts = this->shapeNames.N;
    calcAllDistances();
    filterNeighborMatrixByDist(0.05);
    setGroundParts();
    setNeighborList();//needed?
    setNeccessaryNeighbors(0.05); // indicates if a neighbor's position is lower
  }


  // add a floor to the model
  void addFloor(){
    auto *base = C.addFrame("floor", "");
    base->setShape(rai::ST_ssBox, {18., 18., .3, .02});
    base->setPosition({0.,0.,-.15});
    base->setContact(1.);
  }

  StringA getNeighbours(const rai::String part, const bool floorSet = false) const {
    StringA tmp = {part};
    return getNeighbours(tmp, floorSet);
  }

  // query all Neighbors of a list of parts
  StringA getNeighbours(StringA placed, const bool floor_set=true) const {
    if(floor_set) placed.append("floor");

    uintA placed_idx,nei_idx,ground_idx;
    StringA nei_list;
    //find placed in shapeNames, append idx where found
    for(rai::String s: placed){
      uint idx = 0;
      for(rai::String c: this->shapeNames){
        if(s == c) placed_idx.append(idx);
        idx++;
      }
    }
    if(this->groundParts.N && floor_set){
      ground_idx = getPartIndices(this->groundParts);
      for(uint g:ground_idx){
        if(!placed_idx.contains(g)) nei_idx.append(g);
      }
    }

    for(uint p:placed_idx){
      for(uint i=0;i<this->shapeNames.N;i++){
        //if i is a neighbor of p and i has not been added yet and i has not been placed yet and i is not the floor
        if(!!this->neighborMatrix(p,i) && !nei_idx.contains(i) && !placed_idx.contains(i) && this->shapeNames(i)!="floor"){
          nei_idx.append(i);
        }
      }
    }

    for(uint i:nei_idx){
      nei_list.append(this->shapeNames(i));
    }
    return nei_list;
  }

  // Set ground parts (placeable every round)
  void setGroundParts(){
    StringA tmp = {};
    this->groundParts = this->getNeighbours(tmp);
  }

  // Calculates all pairwise distances in O(N_parts^2)
  void calcAllDistances(){
    C.swiftDelete();
    C.swift();
    rai::Proxy p;
    doubleA distMat(this->shapeNames.N,this->shapeNames.N);
    for(uint i=0;i<this->shapeNames.N;i++){
      for(uint j=i+1;j<this->shapeNames.N;j++){
        p.a = C.getFrame(this->shapeNames(i));
        p.b = C.getFrame(this->shapeNames(j));
        //std::cout<<"try "<<this->shapeNames(i)<< " and "<<this->shapeNames(j)<<std::endl;

        p.calc_coll();
        //std::cout << "dist="<<p.d<<std::endl;

        distMat(i,j)=p.d;
        distMat(j,i)=p.d;
      }
    }
    this->distanceMatrix = distMat;
  }

  // Sets parts as neighbors if sufficiently close
  void filterNeighborMatrixByDist(const double thres=0.05){
    intA neiMat(N_parts,N_parts);
    neiMat = 0;
    for(uint i=0;i<N_parts;i++){
      for(uint j=i+1;j<N_parts;j++){
        if(this->distanceMatrix(i,j)<=thres){
          neiMat(i,j) = 1;
          neiMat(j,i) = 1;
        }
      }
    }
    this->neighborMatrix = neiMat;
  }

  // Initialize list of shape-names in current configuration
  void setShapeNames(){
    StringA names;
    for(rai::Frame* f:C.frames){
      if(f->shape && !f->name.startsWith("_")) names.append(f->name);
    }
    this->shapeNames = names;
  }

  void setNeighborList(){
    intA neiMat(N_parts-1,N_parts-1);
    for(uint i=0;i<N_parts-1;i++){
      for(uint j=0;j<N_parts-1;j++){
        neiMat(i,j) = this->neighborMatrix(i,j);
      }
    }
    uint max_neis=0;
    for(uint i=0;i<N_parts-1;i++){
      uint n_neis = sum(neiMat.row(i));
      if(n_neis>max_neis) max_neis=n_neis;
    }
    intA neiList(N_parts-1,max_neis);
    neiList = -1;
    for(uint i=0;i<N_parts-1;i++){
      for(uint j=0;j<N_parts-1;j++){
        if(neiMat(i,j)){
          uint idx = 0;
          for(int n:neiList.row(i)){
            if(n==-1){
              neiList(i,idx)=j;
              break;
            }
            idx++;
          }
        }
      }
    }
    this->neighborList = neiList;
  }

  // set neccessaryNeighbors
  void setNeccessaryNeighbors(const double tol=0.05){
    intA necMat(this->N_parts,this->N_parts);
    necMat = 0;

    for(uint i=0;i<this->N_parts;i++){
      for(uint j=0;j<this->N_parts;j++){
        if(this->neighborMatrix(i,j)==1){
          rai::Frame* part_i = C.getFrame(this->shapeNames(i));
          rai::Frame* part_j = C.getFrame(this->shapeNames(j));

          const double z_i = part_i->getPosition()(2);
          const double z_j = part_j->getPosition()(2);
          
          if(z_i > z_j && z_i-z_j>tol){
            necMat(i,j) = 2;
          }
          else{
            necMat(i,j) = 1;
          }
        }
      }
    }
    this->neccessaryNeighbors = necMat;
  }

  // delete parts with specific indices, keep the rest
  void deleteSpecific(const uintA &del){
    StringA all_shapes = this->shapeNames;
    StringA del_shapes = {};
    for(uint d:del){
      del_shapes.append(all_shapes(d));
    }
    for(rai::Frame* f: C.getFrames(del_shapes)){
      std::cout << "Frame " << f->name << " will now be deleted!" << std::endl;
      delete f;
    }
    this->setShapeNames();
    this->N_parts = this->shapeNames.N;
  }

  uintA getPartIndices(const StringA &names)const {
    // TODO: replace by find
    uintA indices;
    uint idx = 0;
    for(rai::String s: this->shapeNames){
      if(names.contains(s)) indices.append(idx);
      idx++;
    }
    return indices;
  }

  uint getPartIndex(const rai::String &name)const {
    // TODO: replace by find
    uint idx=0;
    for(const rai::String &s: this->shapeNames){
      if(s==name){
        return idx;
      }
      idx++;
    }
  }
};
