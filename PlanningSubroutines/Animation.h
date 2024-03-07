#pragma once

#include <Kin/kin.h>
#include <Kin/viewer.h>
#include <GL/gl.h>
#include <Gui/opengl.h>

namespace rai {

struct DrawPaths : GLDrawer {
  arr& X;
  DrawPaths(arr& X): X(X) {}
  void glDraw(OpenGL& gl) {
//#ifdef RAI_GL
    glColor(0., 0., 0.);
    for(uint i=0; i<X.d1; i++) {
      glBegin(GL_LINES);
      glLineWidth(2.0);
      for(uint t=0; t<X.d0; t+=2) {
        rai::Transformation pose;
        pose.set(&X(t, i, 0));
//          glTransform(pose);
        glVertex3d(pose.pos.x, pose.pos.y, pose.pos.z);
      }
      glEnd();
    }
//#endif
  }
};

struct Animation{
  struct AnimationPart{
    StringA frameNames;
    uintA frameIDs;
    arr frameCols;
    arr X;
      
    double start = 0.;

    void write(ostream& os) const;
    void read(istream& is);
  };

  rai::Array<AnimationPart> A;

  FrameL prePlannedFrames = {};
  uint tPrePlanned = 0;

  void write(ostream& os) const;
  void read(istream& is);

  uint getT() const;

  mutable std::unordered_map<uint, arr> poseMap;
  void setToTime(rai::Configuration& C, const double t, const double tIgnoreBefore=0) const;
  uintA getActiveFramesAtTime(const double t);

  void play(rai::Configuration& C, bool pause, bool exp=false);
  void drawTrajectories(rai::Configuration& C);
};
stdPipes(Animation::AnimationPart)
stdPipes(Animation)

}

using AnimationPtr = std::shared_ptr<rai::Animation>;
using AnimationPartL = rai::Array<rai::Animation::AnimationPart>;
