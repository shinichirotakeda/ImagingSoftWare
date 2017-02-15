#include "TwoHitComptonEvent.hh"


TwoHitComptonEvent::TwoHitComptonEvent()
{
  m_bCalc = false;
  m_H1Id = 0;
  m_H1DetId = 0;
  m_H1Process = 0;
  m_H1Time = 0.0;
  m_H1PosX = 0.0;
  m_H1PosY = 0.0;
  m_H1PosZ = 0.0;
  m_H1Energy = 0.0;
  m_H2Id = 0;
  m_H2DetId = 0;
  m_H2Process = 0;
  m_H2Time = 0.0;
  m_H2PosX = 0.0;
  m_H2PosY = 0.0;
  m_H2PosZ = 0.0;
  m_H2Energy = 0.0;
}


double TwoHitComptonEvent::CosThetaE() 
{ 
  if (!m_bCalc) { calc(); }
  return m_CosThetaE;
}
  

TVector3 TwoHitComptonEvent::ConeAxis()
{
  if (!m_bCalc) { calc(); }
  return m_ConeAxis;
}


TVector3 TwoHitComptonEvent::ConeVertex()
{
  if (!m_bCalc) { calc(); }
  return m_ConeVertex;
}


void TwoHitComptonEvent::Swap()
{
  m_bCalc = false;

  int H1Id = m_H1Id;
  int H1DetId = m_H1DetId;
  int H1Process = m_H1Process;
  double H1Time = m_H1Time;
  double H1PosX = m_H1PosX;
  double H1PosY = m_H1PosY;
  double H1PosZ = m_H1PosZ;
  double H1Energy = m_H1Energy;
  m_H1Id = m_H2Id;
  m_H1DetId = m_H2DetId;
  m_H1Process = m_H2Process;
  m_H1Time = m_H2Time;
  m_H1PosX = m_H2PosX;
  m_H1PosY = m_H2PosY;
  m_H1PosZ = m_H2PosZ;
  m_H1Energy = m_H2Energy;
  m_H2Id = H1Id;
  m_H2DetId = H1DetId;
  m_H2Process = H1Process;
  m_H2Time = H1Time;
  m_H2PosX = H1PosX;
  m_H2PosY = H1PosY;
  m_H2PosZ = H1PosZ;
  m_H2Energy = H1Energy;
}


void TwoHitComptonEvent::calc()
{
  const double MeV = 1.0;
  const double ELECTRON_MASS = 0.5109989 * MeV;
  m_CosThetaE = 1. - (ELECTRON_MASS*m_H1Energy)/(m_H2Energy*(m_H1Energy+m_H2Energy));
  
  TVector3 v1(m_H1PosX, m_H1PosY, m_H1PosZ);
  TVector3 v2(m_H2PosX, m_H2PosY, m_H2PosZ);
  m_ConeVertex = v1;
  v1 -= v2;
  m_ConeAxis = v1.Unit();
  
  m_bCalc = true;
}

