#ifndef TwoHitComptonEvent_h
#define TwoHitComptonEvent_h

#include <iostream>
#include <vector>

#include "TVector3.h"

class TwoHitComptonEvent
{
public:
  TwoHitComptonEvent();
  virtual ~TwoHitComptonEvent(){}

  int getH1Id(){ return m_H1Id; }
  int getH1DetId() { return m_H1DetId; }
  int getH1Process() { return m_H1Process; }
  double getH1Time() { return m_H1Time; }
  double getH1PosX() { return m_H1PosX; }
  double getH1PosY() { return m_H1PosY; }
  double getH1PosZ() { return m_H1PosZ; }
  double getH1Energy() { return m_H1Energy; }

  int getH2Id(){ return m_H2Id; }
  int getH2DetId() { return m_H2DetId; }
  int getH2Process() { return m_H2Process; }
  double getH2Time() { return m_H2Time; }
  double getH2PosX() { return m_H2PosX; }
  double getH2PosY() { return m_H2PosY; }
  double getH2PosZ() { return m_H2PosZ; }
  double getH2Energy() { return m_H2Energy; }

  void setH1Id(int r){ m_H1Id = r; m_bCalc = false; }
  void setH1DetId(int r){ m_H1DetId = r; m_bCalc = false; }
  void setH1Process(int r){ m_H1Process = r; m_bCalc = false; }
  void setH1Time(double r) { m_H1Time = r; m_bCalc = false; }

  void setH1Pos(double x, double y, double z)
  { 
    m_H1PosX = x;  m_H1PosY = y;  m_H1PosZ = z;
    m_bCalc = false; 
  }

  void setH1Energy(double r) { m_H1Energy = r; m_bCalc = false; }

  void setH2Id(int r){ m_H2Id = r; m_bCalc = false; }
  void setH2DetId(int r){ m_H2DetId = r; m_bCalc = false; }
  void setH2Process(int r){ m_H2Process = r; m_bCalc = false; }
  void setH2Time(double r) { m_H2Time = r; m_bCalc = false; }

  void setH2Pos(double x, double y, double z)
  { 
    m_H2PosX = x;  m_H2PosY = y;  m_H2PosZ = z;
    m_bCalc = false; 
  }

  void setH2Energy(double r) { m_H2Energy = r; m_bCalc = false; }

  double CosThetaE();
  TVector3 ConeAxis();
  TVector3 ConeVertex();
  double TotalEnergy() { return m_H1Energy + m_H2Energy; }
  bool TimeOrder() { return (m_H1Time <= m_H2Time) ? true : false; }
  bool EnergyOrder() { return (m_H1Energy <= m_H2Energy) ? true : false; }
  void Swap();

private:
  bool m_bCalc;

  int m_H1Id;
  int m_H1Process;
  int m_H1DetId;
  double m_H1Time;
  double m_H1PosX;
  double m_H1PosY;
  double m_H1PosZ;
  double m_H1Energy;
  int m_H2Id;
  int m_H2DetId;
  int m_H2Process;
  double m_H2Time;
  double m_H2PosX;
  double m_H2PosY;
  double m_H2PosZ;
  double m_H2Energy;

  double m_CosThetaE;
  TVector3 m_ConeVertex;
  TVector3 m_ConeAxis;

private:
  void calc();
};


#endif // TwoHitComptonEvent_h
