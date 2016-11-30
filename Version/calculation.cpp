#include<iostream>
#include<math.h>
#include<stdio.h>
#include <winsock2.h>
#include<fstream>
#include<strstream>
#include<vector>
#include<stdlib.h>
using namespace std;

//float calculate_Temperature(float pH2,float pCO2,float pCH4,float T,float specific_growth_rate,int Hr)
float calculate_Temperature(double pH2,double pCO2,double pCH4,float T,double specific_growth_rate,double m0,float detaday,float year,float detayear,FILE *ww, FILE *md)
{double pH2_temp, pCO2_temp, pCH4_temp,day_temp;
pH2_temp=pH2;
pCO2_temp=pCO2;
pCH4_temp=pCH4;
int Hr=24;
day_temp=Hr/24;
//����CO2��H2��CO2�ӷ�ѹ���ܽ�ȵ�ת��
//float pH2, pCO2, pCH4;//define partial pressure of H2, CO2, CH4,(unit in atm)
const double Hh2=0.00078,Hch4=0.0006364;//define henry constant of h2 and ch4,(unit in mole/L*atm)
float cH2, cCO2, cCH4;//define concentration of H2, CO2, CH4, (unit in mole/L)
cH2=Hh2*pH2;//calculate the concentration of H2, CO2, CH4 by partial pressure and henry constant
cCO2=pCO2;
cCH4=Hch4*pCH4;
//ͨ��Ũ�ȼ���Q
float Q;
Q=cCH4/cCO2/pow(cH2,4);
//ͨ��Q����detaG1
int detaG0=-131;//define the standard change of gibbs energy(unit in KJ/mole)
const double R=0.008314;//define the gas constant(unit in KJ/mole/K)
float detaG1;//define the change of gibbs energy
//read methanogens-data.txt

//float T;// define the temperature


//float ln(float Q);
float lnQ;
lnQ=log(Q);
detaG1=detaG0+R*T*lnQ;//calculate the change of gibbs energy by formula
//����������������
//int Hr=24;//define time scale 24 hours 

//float specific_growth_rate;//define specific_growth_rate//===doubt===���¶��йأ��Ƿ�õ���һ���ⲿ��.txt�ļ�?
float test=log(2.7);
float td=log(2)/specific_growth_rate;//calculate double growth time
//float day_times;//define the relationship of one day and double growth time
float day_times=Hr/td;//calculate the relationship of one day and double growth time
//float biomass_times;//define how times methanogens biomass increase in one day
float biomass_times=pow(2,day_times);//calculate how times methanogens biomass increase in one day
//===doubt===�о�Ӧ�üӸ�ѭ����������?!
//float m0=0.00001;//define the initial methanogens biomass
float mglobal=2.57e17;//define the global biomass which limits the mb
double mb,mb_temp;
mb=m0*biomass_times; //���������ղ���
mb_temp=mb;


//����K2
double rtns2K;//defiene how many times does reaction of methanogens biomass happen per unit time.(rnts/s)
float Mmole=23.28; // mole mass of methanogens biomass,(unit in g/mole)
float Na=6.02e23;//avogadro constant
rtns2K=(mb-m0)*Na/Mmole/Hr;//calculate how many times does reaction of methanogens biomass happen per unit time.(rnts/s)

//ͨ������ƽ�����K1
float rtns1K;//defiene how many times does reaction of methane production happen per unit time.(rnts/s)
int detaG2=-50;// the energy for cell production(unit in KJ/mole)  
rtns1K=rtns2K*detaG2/detaG1;//calculate how many times does reaction of methane production happen per unit time.(rnts/s)

//ͨ��K1����mch4�ղ�����nch4,����Vch4��Ȼ�����pCH4
float nCH4,Vmole,VCH4;
nCH4=rtns1K*Hr/Na;
Vmole=T/273*22.4/3;
VCH4=nCH4*Vmole;

//����VH2��VCO2
float nb,nH2,nCO2,VH2,VCO2;
nb=mb/Mmole;
nH2=2.09*nb+4*nCH4;
nCO2=nb+nCH4;
VH2=nH2*Vmole;
VCO2=nCO2*Vmole;

//��������е�pCH4(�����������CH4�����ۻ�)��pH2��pCO2��������ԭ�еļ�ȥ��������ɵģ�
float Va=4e21;
float pCH4_new, pH2_new,pCO2_new;
pCH4_new=VCH4/Va;
pH2_new=VH2/Va;
pCO2_new=VH2/Va;
pCH4=pCH4_new+pCH4_temp;
pH2=pH2_temp-pH2_new;
pCO2=pCO2_temp-pCO2_new;

float i,Tdegree, Carbon, Hydrogen;
ww=fopen("C:\\methanogens_atmosphere_model\\a.txt","r");

double minlenght,minlenght_temp; 
double squrepCO2,squrepH2;

squrepCO2=pow(0.031623,2);
squrepH2=pow(0.1,2);
minlenght_temp=squrepCO2+squrepH2;


double pH2_temp_min_length=0, pCO2_temp_min_length=0;

//while(!feof(ww)){
while(1){
//while(fscanf(ww,"%f %f %f",&Tdegree,&Carbon,&Hydrogen)==2){
fscanf(ww,"%f %f %f %f", &i, &Tdegree,&Hydrogen,&Carbon);
	if (feof(ww)) 
		break;
	   squrepCO2=pow((pCO2-Carbon),2);
       squrepH2=pow((pH2-Hydrogen),2);
	   minlenght=squrepCO2+squrepH2;
	   //minlenght_temp=minlenght;
		if(minlenght < 0)
			return T;
	   if(minlenght<minlenght_temp)
	   {minlenght_temp=minlenght;
		   //pCO2=Carbon;
	   pH2_temp_min_length = Hydrogen;
	   //pH2=Hydrogen;
	   pCO2_temp_min_length =  Carbon;
	   T=Tdegree;} 


	   	if(pH2_temp_min_length==0)
			return T;

	   //printf(" %f %f %f %f\n",  Carbon, Hydrogen, T, minlenght);
	   //else{return pCO2, pH2, T;}
}


//intf("The right point:  %f %f %f %f\n",  pCO2_temp_min_length, pH2_temp_min_length, T, minlenght);
	   //pH2=pH2_temp_min_length;
	    
	  // pCO2=pCO2_temp_min_length;

	   FILE *right_point_file;
right_point_file=fopen("C:\\methanogens_atmosphere_model\\right_point_file.txt","a");
//fprintf(right_point_file, "The right point:  %f %f %f %f\n",  pCO2_temp_min_length, pH2_temp_min_length, T, minlenght);
	   fclose(right_point_file); 

fclose(ww);

//float mname;
float temperature, growthrate;
md=fopen("C:\\methanogens_atmosphere_model\\a_b.txt","r");
while(!feof(md)){


	fscanf(md,"%f %f\n",&temperature,&growthrate);

	if (feof(md)) break;

//printf("%f %f\n", temperature, growthrate);
if(T==temperature)
specific_growth_rate=growthrate;//Ӧ�ö����Ƕ�Ӧ�е�growthrate
}
fclose(md);
if(mb<=mglobal)
{float m_temp;
	m_temp=mb;
	m0=m_temp;
	//return mb;
	//detaday++;
	  }
      else{
	  //mb=m0;
		  m0 = 1e-5;
	
		  detayear++;
}

m0=mb_temp;
float day;
day=day_temp+detaday;
year=day/365+detayear;
FILE *rt;
rt=fopen("C:\\methanogens_atmosphere_model\\result.txt","a");
fprintf(rt,"%f %f %f %f %f %f %f\n",detayear,T,pH2,pCO2,pCH4,mb,m0);
printf("%f %f %f %f %f %f %f\n",detayear,T,pH2,pCO2,pCH4,mb,m0);
fclose(rt);
if(T<=277)
{
return T;
}else
{
calculate_Temperature(pH2,pCO2,pCH4,T,specific_growth_rate,m0,detaday,year,detayear,ww,md);
}
//double day; //define time, hr=one hour, day=24hrs, year=365 days

//printf("%f %f %f  \n",year,day,Hr);
//day+=Hr/24;

//���Ӵ�ӡ���֣�print��
//T=(��wordsworth��data)


//pCH4+=pCH4;
//pH2-=pH2;
//pCO2-=pCO2;
//FILE *ww, *md, *a, *b, *a_b;

/*if((md=fopen("C:\methanogens_atmosphere_model\methanogen-data.txt","w"))==NULL)
{printf("open file fail\n");
return;
}*/
}
void main()
{
//int n,i;

float pH2=0.1;
float pCO2=0.031623;
float pCH4=1e-6;
float T=285;
double specific_growth_rate=0.09354213;
//int Hr=24;
double m0=1e-5;
//float day=0;
float detaday=1;
float year=0;
float detayear=0;
FILE *ww, *md;
//float mname;


calculate_Temperature(pH2,pCO2,pCH4,T,specific_growth_rate,m0,detaday,year,detayear,ww,md);

//calculate_Temperature(pH2,pCO2,pCH4, FILE *ww,FILE *md);===?===��α�ʾ���ŵ��ú�����������
} 

