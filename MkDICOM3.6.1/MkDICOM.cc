#include "MkDICOM.hh"

#include "NextCLI.hh"


#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include "time.h"

using namespace anl;

void MkDICOM::Init(){
  
  std::string name1 = "input.root";
  CLread("input root file name ",&name1);
  TFile *fin = new TFile(name1.c_str());
  std::string name11 = "image_000";
  CLread("hist name",&name11);
  TH3D *hin = (TH3D *)fin->Get(name11.c_str());
  
  Set3DHist(hin);

  name_folder = "COMPTON";
  CLread("Folder name ",&name_folder);

  
  std::string name2 ="TAKEDA";
  CLread("Patient Name (6 words) ",&name2);
  p_name = name2;
  //  std::cout << p_name << std::endl;
  
  time_t now = time(NULL);
  struct tm *pnow = localtime(&now);
  std::ostringstream year;
  year << std::setfill('0') << std::setw(6) << 1900+pnow->tm_year;
  std::ostringstream mon;
  mon << std::setfill('0') << std::setw(2) << pnow->tm_mon;
  std::ostringstream day;
  day << std::setfill('0') << std::setw(2) << pnow->tm_mday;
  std::ostringstream hour;
  hour << std::setfill('0') << std::setw(2)<< pnow->tm_hour;
  std::ostringstream min;
  min << std::setfill('0') << std::setw(2) << pnow->tm_min;
  std::ostringstream sec;
  sec << std::setfill('0') << std::setw(2) << pnow->tm_sec;  
  //  std::cout << year.str() << mon.str() << day.str() << hour.str() << min.str() << sec.str() << std::endl;

  std::string name4 = sec.str() + min.str() + (hour.str()).at(0) + "."
    + (hour.str()).at(1) + day.str() + (mon.str()).at(0) + "." + (mon.str()).at(1) + year.str();
    
  //  CLread("Study Number ",&name4);
  n_study = name4;
  
  std::string name3 = sec.str() + min.str() + (hour.str()).at(0) + "." + (hour.str()).at(1) + "." + (day.str()).at(0);
  //  CLread("SOP Instance Number ",&name3);
  base_name = name3;
  
  //  std::string sli="1";
  //  CLread("Slice ? [mm]",&sli);
  //  slice = sli;
  
  sli_dir = 1;
  CLread("Slice Direction ? X(0), Y(1), Z(2)", &sli_dir);
  
  namespace fs = boost::filesystem;
  fs::path dir(name_folder.c_str());
  fs::create_directory(dir);
  
}


void MkDICOM::Set3DHist(TH3 *h3)
{

  sbp= h3;
  double max,min;
  
  max = 0.0;
  min = 10000.0;
  
  for(int i=1;i<=sbp->GetXaxis()->GetNbins();i++){
    for(int j=1;j<=sbp->GetYaxis()->GetNbins();j++){
      for(int k=1;k<=sbp->GetZaxis()->GetNbins();k++){
        int gbin = sbp->GetBin(i,j,k);
        double v = sbp->GetBinContent(gbin);
	
        if(max<v){
          max = v;
        }
	
        if(v<min){
          min = v;
        }
      }
    }
  }
    
  h_maxvalue = max;
  h_minvalue = min;
  
  
}

void MkDICOM::DICOMFileGen()
{

  int nfile = ChartStackVec.size();
  for(int i=0;i<nfile;i++){
    DICOMFileFactory(ChartStackVec[i]);
  }
  
}


void MkDICOM::DICOMFileFactory(ChartStack *s)
{
  
  std::string base = "./" + name_folder + "/";
  std::string name(s->Get2DHist()->GetName());
  std::string filename= base+name;
  filename = filename +".dcm";
  std::cout << filename << std::endl;

  char c_name[256];
  strcpy(c_name,filename.c_str());
  
  std::ofstream fout;
  fout.open(c_name,std::ios::binary);
  if(!fout){
    std::cout << "Can not open file " <<std::endl;
  }
  
  
  //  std::cout << "sizeof(unsigned short int) " << sizeof(unsigned short int) << std::endl;
  //  std::cout << "sizeof(unsigned int) " << sizeof(unsigned int) << std::endl;
  //  std::cout << "sizeof(int) " << sizeof(int) << std::endl;
  //  std::cout << "sizeof(long) " << sizeof(long) << std::endl;
  //  std::cout << "sizeof(char) " << sizeof(char) << std::endl;

  AddZeroHedder(&fout,s);
  Add0002Tag(&fout,s);
  Add0008Tag(&fout,s);
  Add0010Tag(&fout,s);
  Add0018Tag(&fout,s);
  Add0020Tag(&fout,s);
  Add0028Tag(&fout,s);
  Add0032Tag(&fout,s);
  Add0040Tag(&fout,s);
  AddPixelData(&fout,s);
  
  fout.close();

}


void MkDICOM::AddZeroHedder(ofstream *f,const ChartStack *s)
{
  //128 byte Null
  char a[128];
  for(int i=0;i<128;i++){
    a[i]=0x00;
  }
  f->write(a,128);
}

void MkDICOM::Add0002Tag(ofstream *f,const ChartStack *s)
{
  
  unsigned short int group = 0x0002;
  unsigned short int element;
  
  unsigned short int data16; // 16bit size
  unsigned int data32; // 32bit size
  
  char dicom[]="DICM";
  char ul[]="UL";
  char ob[]="OB";
  char sh[]="SH";
  char ae[]="AE";
  char ui[]="UI";
  
  f->write(dicom,4);

  //0002 0000 
  element = 0x0000;
  data16 = 0x0004;
  data32 = 0x000000c0; // 192 ??
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ul,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  f->write((char *)&data32, sizeof(unsigned int));
  
  //0002 0001
  element = 0x0001;
  data16 = 0x0000;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ob,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 = 0x0002;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16=0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16=0x0100;
  f->write((char *)&data16, sizeof(unsigned short int));

  //0002 0002
  element = 0x0002;
  data16 = 0x001a;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ui,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  char mssopcuid[]="1.2.840.10008.5.1.4.1.1.2";
  f->write(mssopcuid,data16);

  //0002 0003
  element = 0x0003;
  data16 = 0x0030;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ui,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  std::string t_base="1.3.12.2.1107.5.1.4.";
  t_base = t_base + base_name;
  std::string str = s->Get2DHist()->GetName();
  std::string s_str = str.substr(7,3);
  t_base = t_base + ".000000000000000";
  t_base = t_base + s_str;
  char mssopiuid[t_base.size()];
  strcpy(mssopiuid,t_base.c_str());
  f->write(mssopiuid,data16);

  //0002 0010
  element = 0x0010;
  data16 = 0x0012;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ui,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  char tsuid[]="1.2.840.10008.1.2";
  f->write(tsuid,data16);

  //0002 0012
  element = 0x0012;
  data16 = 0x0016;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ui,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  char icuid[]="1.3.6.1.4.1.19291.2.1";
  f->write(icuid,data16);

  //0002 0013
  element = 0x0013;
  data16 = 0x000a;
  char osirix001[] = "OSIRIX001";
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(sh,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  f->write(osirix001,10);

  //0002 0016
  element = 0x0016;
  data16 = 0x0006;
  char osirix[] = "OSIRIX";
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write(ae,2);
  f->write((char *)&data16, sizeof(unsigned short int));
  f->write(osirix,6);

}

void MkDICOM::Add0008Tag(ofstream *f,const ChartStack *s)
{


  unsigned short int group = 0x0008;
  unsigned short int element;
  unsigned int data32; // 32bit size
  unsigned short int data16; // 16bit size

  std::string id_study = "1.2.840.113745.101000.1008000." + n_study;


  //0008 0005
  element = 0x0005;
  data32 = 0x0000000a;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char iso_ir[]="ISO_IR 100";
  f->write(iso_ir,10);

  //0008 0008
  element = 0x0008;
  data32 = 0x00000034; // 52 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char imagetype[]="ORIGINAL\\PRIMARY\\AXIAL\\CT_SOMS SPI\\HEADER_CORRECTED";
  f->write(imagetype,52);

  //0008 0016
  element = 0x0016;
  data32 = 0x0000001a; // 26 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char sopclassuid[]="1.2.840.10008.5.1.4.1.1.2";
  f->write(sopclassuid,26);

  //0008 0018


  element = 0x0018;
  data32 = 0x00000030; // 48 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  std::string t_base="1.3.12.2.1107.5.1.4.";
  t_base = t_base + base_name;
  std::string str = s->Get2DHist()->GetName();
  std::string s_str = str.substr(7,3);
  t_base = t_base + ".000000000000000";
  t_base = t_base + s_str;
  char sopinstanceuid[t_base.size()];
  strcpy(sopinstanceuid,t_base.c_str());
  f->write(sopinstanceuid,48);

  time_t now = time(NULL);
  struct tm *pnow = localtime(&now);
  std::ostringstream year;
  year << std::setfill('0') << std::setw(4) << 1900+pnow->tm_year;
  std::ostringstream mon;
  mon << std::setfill('0') << std::setw(2) << pnow->tm_mon;
  std::ostringstream day;
  day << std::setfill('0') << std::setw(2) << pnow->tm_mday;
  std::ostringstream hour;
  hour << std::setfill('0') << std::setw(2)<< pnow->tm_hour;
  std::ostringstream min;
  min << std::setfill('0') << std::setw(2) << pnow->tm_min;
  std::ostringstream sec;
  sec << std::setfill('0') << std::setw(2) << pnow->tm_sec;  
  
  //0008 0020
  element = 0x0020;
  data32 = 0x00000008; // 8 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  //  char date[] ="20000101";
  std::string date = year.str() + mon.str()+ day.str();
  f->write(date.c_str(),8);

  //0008 0021
  element = 0x0021;
  data32 = 0x00000008; // 8 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(date.c_str(),8);

  //0008 0022
  element = 0x0022;
  data32 = 0x00000008; // 8 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(date.c_str(),8);

  //0008 0023
  element = 0x0023;
  data32 = 0x00000008; // 8 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(date.c_str(),8);

  //0008 0030
  element = 0x0030;
  data32 = 0x0000000e; // 14 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  //  char time[]="100010.200000 ";
  std::string time = hour.str() + min.str()+ sec.str() + ".200000 ";

  f->write(time.c_str(),14);

  //0008 0031
  element = 0x0031;
  data32 = 0x0000000e; // 14 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(time.c_str(),14);

  //0008 0032
  element = 0x0032;
  data32 = 0x0000000e; // 14 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(time.c_str(),14);

  //0008 0033
  element = 0x0033;
  data32 = 0x0000000e; // 14 byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(time.c_str(),14);

  //0008 0050
  element = 0x0050;
  data32 = 0x00000008; // 8byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char AccessNum[]="0000001";
  f->write(AccessNum,8);

  //0008 0060
  element = 0x0060;
  data32 = 0x00000002; // 2byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char mod[]="CC";
  f->write(mod,2);

  //0008 0070
  element = 0x0070;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char manu[]="RIKEN";
  f->write(manu,6);

  //0008 0080
  element = 0x0080;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(manu,6);

  //0008 0081
  element = 0x0081;
  data32 = 0x00000012; // 18byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char address[]="TYUOKU/KOBE/HYOGO";
  f->write(address,18);

  //0008 0090
  element = 0x0090;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char phname[]="TAKEDA";
  f->write(phname,6);

  //0008 1010
  element = 0x1010;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char station[]="CC0001";
  f->write(station,6);

  //0008 1030
  element = 0x1030;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char study[]="CCTEST";
  f->write(study,6);

  //0008 103e
  element = 0x103e;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char series[]="CCTEST";
  f->write(series,6);

  //0008 1090
  element = 0x1090;
  data32 = 0x00000006; // 6byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char model[]="CCTEST";
  f->write(model,6);

  //0008 1110
  element = 0x1110;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1150;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000018;// 24 byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test1[]="1.2.840.10008.3.1.2.3.1";
  f->write(test1,24);
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1155;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000030;// 48 byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test2[50];
  strcpy(test2,id_study.c_str());
  f->write(test2,48);
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe00d;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe0dd;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));

  //0008 1111
  element = 0x1111;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1150;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000018;// 24 byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test3[]="1.2.840.10008.3.1.2.3.2";
  f->write(test3,24);
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1155;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000030;// 48 byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test4[50];
  strcpy(test4,id_study.c_str());
  f->write(test4,48);
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe00d;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe0dd;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));

  //0008 1120
  element = 0x1120;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1150;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000018;// 24 byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test5[]="1.2.840.10008.3.1.2.1.1";
  f->write(test5,24);
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1155;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000030;// 48 byte
  f->write((char *)&data32, sizeof(unsigned int));
  //  char test6[]="1.2.840.113745.101000.1008000.37317.4904.3609243";
  char test6[50];
  strcpy(test6,id_study.c_str());
  f->write(test6,48);
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe00d;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe0dd;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));

  //0008 1140
  element = 0x1140;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1150;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x0000001a;// 26 byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test7[]="1.2.840.10008.5.1.4.1.1.2";
  f->write(test7,26);
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1155;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000030;// 48 byte
  f->write((char *)&data32, sizeof(unsigned int));
  //  char test8[]="1.3.12.2.1107.5.1.4.36085.4.0.13454957584835538";
  std::string tt_base="1.3.12.2.1107.5.1.4.";
  tt_base = tt_base + base_name;
  tt_base = tt_base + ".13454957584835538";
  char test8[48];
  strcpy(test8,tt_base.c_str());
  f->write(test8,48);
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe00d;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe0dd;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));

  //0008 2112
  element = 0x2112;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xffff;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1150;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000014;//  20byte
  f->write((char *)&data32, sizeof(unsigned int));
  char test9[]="1.3.12.2.1107.5.9.1";
  f->write(test9,20);
  data16 =0x0008;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x1155;
  f->write((char *)&data16, sizeof(unsigned short int));
  data32 =0x00000030;// 48 byte
  f->write((char *)&data32, sizeof(unsigned int));
  //  char test10[]="1.3.12.2.1107.5.1.4.36085.4.0.13457317828212953";
  std::string ttt_base="1.3.12.2.1107.5.1.4.";
  ttt_base = ttt_base + base_name;
  ttt_base = ttt_base + ".13457317828212953";
  char test10[48];
  strcpy(test10,ttt_base.c_str());
  f->write(test10,48);
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe00d;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xfffe;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0xe0dd;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
  data16 =0x0000;
  f->write((char *)&data16, sizeof(unsigned short int));
}

void MkDICOM::Add0010Tag(ofstream *f,const ChartStack *s)
{

  unsigned short int group = 0x0010;
  unsigned short int element;
  unsigned int data32; // 32bit size

  //0010 0010
  element = 0x0010; 
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  
  data32 = p_name.size();
  if(data32%2){
    data32++;
  }
  char patiname[p_name.size()];
  strcpy(patiname,p_name.c_str());
  
  //  char patiname[]="TAKEDA";
  //  data32=0x00000006;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(patiname,6);


  //0010 0020
  element = 0x0020;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  
  data32 = p_name.size();
  if(data32%2){
    data32++;
  }
  char patiid[p_name.size()];
  strcpy(patiid,p_name.c_str());
  
  // char patiid[]="TAKEDA";
  //  data32=0x00000006;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(patiid,6);

  //0010 1010
  element = 0x1010;
  data32 = 0x00000004; // 4byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char age[]="000Y";
  f->write(age,4);

  //0010 1030
  element = 0x1030;
  data32 = 0x0000002; // 4byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char weight[]="0";
  f->write(weight,2);

}

void MkDICOM::Add0018Tag(ofstream *f,const ChartStack *s)
{
  unsigned short int group = 0x0018;
  unsigned short int element;
  unsigned int data32; // 32bit size

  //0018 0015
  element = 0x0015;
  data32 = 0x00000004; // 4byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char bpart[]="NECK";
  f->write(bpart,4);

  //0018 0050
  element = 0x0050;
  data32 = 0x00000004; // 2byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  //char slice[]="0.5";
  char sss[4];
  strcpy(sss,slice.c_str());
  f->write(sss,4);

  //0018 0060
  //0018 1020
  //0018 1030
  //0018 1100
  //0018 1110
  //0018 1111
  //0018 1120
  //0018 1130
  //0018 1140
  //0018 1150
  //0018 1151
  //0018 1152
  //0018 1160
  //0018 1170
  //0018 1190
  //0018 1200
  //0018 1201
  //0018 1210
  //0018 5100

}

void MkDICOM::Add0020Tag(ofstream *f,const ChartStack *s)
{
  unsigned short int group = 0x0020;
  unsigned short int element;
  unsigned int data32; // 32bit size
  std::string id_study = "1.2.840.113745.101000.1008000." + n_study;

 
  //0020 000d
  
  element = 0x000d;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data32 =0x00000030;//  48byte
  f->write((char *)&data32, sizeof(unsigned int));  
  char test1[50];
  strcpy(test1,id_study.c_str());
  f->write(test1,48);
  


 //0020 000e
  element = 0x000e;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data32 =0x00000030;//  48byte
  f->write((char *)&data32, sizeof(unsigned int));
  //  char test2[]="1.3.12.2.1107.5.1.4.36085.4.0.13457320123409917";
  std::string t_base="1.3.12.2.1107.5.1.4.";
  t_base = t_base + base_name;
  t_base = t_base + ".13457320123409917";
  char test2[48];
  strcpy(test2,t_base.c_str());
  f->write(test2,48);
  


  //0020 0010
  element = 0x0010;
  data32 = 0x00000008; // 8byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char studyid[]="0000001";
  f->write(studyid,8);

  //0020 0011
  element = 0x0011;
  data32 = 0x00000002; // 2byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char seriesnum[]="1";
  f->write(seriesnum,2);

  //0020 0012
  element = 0x0012;
  data32 = 0x00000002; // 2byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char anum[]="1";
  f->write(anum,2);


  //0020 0013
  element = 0x0013;
  data32 = 0x00000004; // 2byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  std::string str = s->Get2DHist()->GetName();
  std::string s_str = str.substr(7,3);
  std::string ss_str;

  /*
  if(s_str.find("0") == 0){
    ss_str = s_str.substr(1,1); 
  }else{
    ss_str =s_str;
  }
  */
  ss_str = s_str;
  char inum[data32];
  strcpy(inum,ss_str.c_str());
  f->write(inum,data32);


  //0020 0032
  element = 0x0032;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  double xcenter = s->Get2DHist()->GetXaxis()->GetBinCenter(1);
  double ycenter = s->Get2DHist()->GetYaxis()->GetBinCenter(1);
  double slicepos = s->GetSlicePos();
  char x[256];
  char y[256];
  char z[256];
  char buff[256];
  sprintf(x,"%.2f",xcenter);
  sprintf(y,"%.2f",ycenter);
  sprintf(z,"%.2f",slicepos);
  strcpy(buff,x);
  strcat(buff,"\\");
  strcat(buff,y);
  strcat(buff,"\\");
  strcat(buff,z);
  std::string aaa = buff;
  data32 = aaa.size();
  if(data32%2){
    data32++;
  }
  char imagepos[aaa.size()];
  strcpy(imagepos,aaa.c_str());
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(imagepos,data32);




  //0020 0037
  element = 0x0037;
  data32 = 0x0000000c; // 12byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char imageorientation[]="1\\0\\0\\0\\1\\0";
  f->write(imageorientation,12);



  //0020 0052
  element = 0x0052;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  data32 =0x00000030;//  48byte
  f->write((char *)&data32, sizeof(unsigned int));
  // char test3[]="1.3.12.2.1107.5.1.4.36085.4.0.134549532215052935";
  std::string tt_base="1.3.12.2.1107.5.1.4.";
  tt_base = tt_base + base_name;
  tt_base = tt_base + ".134549532215052935";
  char test3[48];
  strcpy(test3,tt_base.c_str());
  f->write(test3,48);


  //0020 1041
  element = 0x1041;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  std::string bbb = z;
  data32 = bbb.size();
  if(data32%2){
    data32++;
  }


  char slicelocation[bbb.size()];
  strcpy(slicelocation,bbb.c_str());
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(slicelocation,data32);


}

void MkDICOM::Add0028Tag(ofstream *f,const ChartStack *s)
{
  unsigned short int group = 0x0028;
  unsigned short int element;

  unsigned short int data16; // 16bit size
  unsigned int data32; // 32bit size




  // 0028 0002
  element = 0x0002;
  data32 = 0x00000002; // 2byte
  data16 = 0x0001;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));
    
  //0028 0004
  element = 0x0004;
  data32 = 0x0000000c; // 12byte
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  char photometric[]="MONOCHROME2";
  f->write(photometric,12);

  //0028 0010
  element = 0x0010;
  data32 = 0x00000002; // 2byte
  data16 = s->Get2DHist()->GetXaxis()->GetNbins();  
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));

  //0028 0011
  element = 0x0011;
  data32 = 0x00000002; // 2byte
  data16 = s->Get2DHist()->GetYaxis()->GetNbins();  
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));

  //0028 0030
  element = 0x0030;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  double x_pwidth = s->Get2DHist()->GetXaxis()->GetBinWidth(1);
  double y_pwidth = s->Get2DHist()->GetYaxis()->GetBinWidth(1);
  char x[256];
  char y[256];
  char buff[256];
  sprintf(x,"%.2f",x_pwidth);
  sprintf(y,"%.2f",y_pwidth);
  strcpy(buff,x);
  strcat(buff,"\\");
  strcat(buff,y);
  std::string a = buff;
  data32 = a.size();
  if(data32%2){
    data32++;
  }
  char pixelspace[a.size()];
  strcpy(pixelspace,a.c_str());
  f->write((char *)&data32, sizeof(unsigned int));
  f->write(pixelspace,data32);

  //0028 0100
  element = 0x0100;
  data32 = 0x00000002; // 2byte
  data16 = 16;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));

  //0028 0101
  element = 0x0101;
  data32 = 0x00000002; // 2byte
  data16 = 16;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));

  //0028 0102
  element = 0x0102;
  data32 = 0x00000002; // 2byte
  data16 = 16;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));

  //0028 0103
  element = 0x0103;
  data32 = 0x00000002; // 2byte
  data16 = 0;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&data16, sizeof(unsigned short int));

  //0028 1050
  element = 0x1050;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  char w_center[]="32768";
  data32=6;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&w_center, data32);

  //0028 1051
  element = 0x1051;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  char w_width[]="65536";
  data32=6;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&w_width, data32);

  //0028 1052
  element = 0x1052;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  char rescale[]="0";
  data32=2;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&rescale, data32);

  //0028 1053
  element = 0x1053;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );

  char slope[]="1";
  data32=2;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&slope, data32);

  //0028 1055
  element = 0x1055;
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  char wcwe[]="WINDOW1\\WINDOW2";
  data32=16;
  f->write((char *)&data32, sizeof(unsigned int));
  f->write((char *)&wcwe, data32);


}

void MkDICOM::Add0032Tag(ofstream *f,const ChartStack *s)
{
 
}

void MkDICOM::Add0040Tag(ofstream *f,const ChartStack *s)
{


}

void MkDICOM::AddPixelData(ofstream *f,const ChartStack *s)
{
  unsigned short int group;
  unsigned short int element;
  unsigned long int datalength;

  group = 0x7FE0;
  element = 0x0010;
  datalength = sizeof(unsigned short int)
    * (s->Get2DHist()->GetXaxis()->GetNbins())
    * (s->Get2DHist()->GetYaxis()->GetNbins());
  
  f->write((char *)&group, sizeof(unsigned short int) );
  f->write((char *)&element, sizeof(unsigned short int) );
  f->write((char *)&datalength, sizeof(unsigned int) );

  unsigned  short int bb;


  for(int i=1;i<=s->Get2DHist()->GetXaxis()->GetNbins();i++){
    for(int j=1;j<=s->Get2DHist()->GetYaxis()->GetNbins();j++){
      double data = s->Get2DHist()->GetBinContent(i,j);
      bb = (unsigned short int) ((double)0xffff * (data-h_minvalue)/(h_maxvalue-h_minvalue));
      f->write((char *)&bb, sizeof(unsigned short int));
    }
  }



}


void MkDICOM::MkChart()
{


  if(sli_dir==1){
    const int nbin = sbp->GetYaxis()->GetNbins();
    const double width = sbp->GetYaxis()->GetBinWidth(1);
    slice = std::to_string(width);
    std::cout << "Thickness " << slice << std::endl;
    const int xmin = sbp->GetXaxis()->GetXmin();
    const int xmax = sbp->GetXaxis()->GetXmax();
    const int xbin = sbp->GetXaxis()->GetNbins();
    
    const int zmin = sbp->GetZaxis()->GetXmin();
    const int zmax = sbp->GetZaxis()->GetXmax();
    const int zbin = sbp->GetZaxis()->GetNbins();
    
    TH2D *hh[nbin];
    char name[256];
    ChartStack *chart[nbin];
    
    for(int i=1;i<=nbin;i++){
      
      chart[i-1] = new ChartStack();
      
      double slicepostion;
      sprintf(name,"COMPTON%03d",i);
      hh[i-1] = new TH2D(name,name,xbin,xmin,xmax,zbin,zmin,zmax);
      slicepostion = sbp->GetYaxis()->GetBinCenter(i);
      
      
      for(int j=1;j<=xbin;j++){
	for(int k=1;k<=zbin;k++){
	  int gbin = sbp->GetBin(j,i,k);
	  hh[i-1]->SetBinContent(j,k,sbp->GetBinContent(gbin));
	}
	
	
      }
      
      chart[i-1]->Set2DHist(hh[i-1]);
      chart[i-1]->SetSliceLocation(slicepostion);
      
      ChartStackVec.push_back(chart[i-1]);
      
    }
  }else if(sli_dir==0){

    const int nbin = sbp->GetXaxis()->GetNbins();
    const double width = sbp->GetXaxis()->GetBinWidth(1);
    slice = std::to_string(width);
    std::cout << "Thickness " << slice << std::endl;
    const int ymin = sbp->GetYaxis()->GetXmin();
    const int ymax = sbp->GetYaxis()->GetXmax();
    const int ybin = sbp->GetYaxis()->GetNbins();
    
    const int zmin = sbp->GetZaxis()->GetXmin();
    const int zmax = sbp->GetZaxis()->GetXmax();
    const int zbin = sbp->GetZaxis()->GetNbins();
    
    TH2D *hh[nbin];
    char name[256];
    ChartStack *chart[nbin];
    
    for(int i=1;i<=nbin;i++){
      
      chart[i-1] = new ChartStack();
      
      double slicepostion;
      sprintf(name,"COMPTON%03d",i);
      hh[i-1] = new TH2D(name,name,ybin,ymin,ymax,zbin,zmin,zmax);
      slicepostion = sbp->GetXaxis()->GetBinCenter(i);
      
      
      for(int j=1;j<=ybin;j++){
	for(int k=1;k<=zbin;k++){
	  int gbin = sbp->GetBin(i,j,k);
	  hh[i-1]->SetBinContent(j,k,sbp->GetBinContent(gbin));
	}
	
	
      }
      
      chart[i-1]->Set2DHist(hh[i-1]);
      chart[i-1]->SetSliceLocation(slicepostion);
      
      ChartStackVec.push_back(chart[i-1]);
      
    }

  }else if(sli_dir==2){

    const int nbin = sbp->GetZaxis()->GetNbins();
    const double width = sbp->GetZaxis()->GetBinWidth(1);
    slice = std::to_string(width);
    std::cout << "Thickness " << slice << std::endl;
    const int ymin = sbp->GetYaxis()->GetXmin();
    const int ymax = sbp->GetYaxis()->GetXmax();
    const int ybin = sbp->GetYaxis()->GetNbins();
    
    const int xmin = sbp->GetXaxis()->GetXmin();
    const int xmax = sbp->GetXaxis()->GetXmax();
    const int xbin = sbp->GetXaxis()->GetNbins();
    
    TH2D *hh[nbin];
    char name[256];
    ChartStack *chart[nbin];
    
    for(int i=1;i<=nbin;i++){
      
      chart[i-1] = new ChartStack();
      
      double slicepostion;
      sprintf(name,"COMPTON%03d",i);
      hh[i-1] = new TH2D(name,name,xbin,xmin,xmax,ybin,ymin,ymax);
      slicepostion = sbp->GetZaxis()->GetBinCenter(i);
      
      
      for(int j=1;j<=ybin;j++){
	for(int k=1;k<=xbin;k++){
	  int gbin = sbp->GetBin(k,j,i);
	  hh[i-1]->SetBinContent(k,j,sbp->GetBinContent(gbin));
	}
	
	
      }
      
      chart[i-1]->Set2DHist(hh[i-1]);
      chart[i-1]->SetSliceLocation(slicepostion);
      
      ChartStackVec.push_back(chart[i-1]);
      
    }


  }else{
    std::cout << "Direction is only 0,1,2" << std::endl;
  } 
    



}


void MkDICOM::PrintChart(){

 TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,1000);
 c1->Divide(8,8);

 for(int i=1;i<=64;i++){

   c1->cd(i);
   ChartStackVec[i-1]->Print2DHist();
   c1->Update();

 }





}





      
      


