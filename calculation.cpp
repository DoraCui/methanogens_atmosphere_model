#include<iostream>
#include<math.h>
#include<stdio.h>
#include <winsock2.h>
#include<fstream>
#include<strstream>
#include<vector>
#include<stdlib.h>
#include <string.h> 
#include <stdio.h>

using namespace std;

void calculate_Temperature(double &pH2,
	double &pCO2,
	double &pCH4,
	float  &T,
	float &specific_growth_rate,
	double &m0,
	int  &detaday,
	double  &year,
	int  &detayear,
	double &rest_m_N,
	int times,
	int runtimes
	)
{
	/*H2每年的新增量*/
	//double H2_vol = 9.25e10;//unit in hydrogen molecules/cm2/s
	double H2_vol = 9.25e10;
	const double Na = 6.02e23, m_at = 5.3e21, S_at = 4 * 3.14*pow(6.4e8, 2), t_year = 3600 * 24 * 365;//avogadro constant
	double ftot, H2_esc, H2_accum, pH2_net;
	
	//carbon source
	double Fmeta;//metamorphic flux=Fhydro/0.9
	const double Fridge = 2e12; //5e12;//5e14; //3.01e14*//2e18;//outgassing at the midocean ridge unit in mol yr - 1
	const double Farc = 2.5e12;//outgassing at arc volcanoed-- mol yr - 1
	//double  Fsio2w, CO2_accumulate, pCO2_net;
	const double Fhydro = 1.65e12;	//carbonatization of oceanic crust in warm hydrothermal systems--unit in mol/yr
	//const double Fej = 3e14;	//carbonatization of impact ejecta--unit in moles yr - 1

	/*natural nitrogen fixation*/
	double NO_energy_yield = 1.3e16;//NO energy yield, unit in moleculars/J
	double light_energy = 1.1e18;//lightning energy yield, unit in J/year
	double initial_NO = 3e11*30/14;//3e11 g N whenever CO2 was the dominant atmospheric gas, *30/14 单位是 g / yr
	//double amount_ini_NO = initial_NO * 7e9;//unit in g 从45亿到38亿年累积起来的NO
	double NO_energy_yield_net = 2e14;//NO energy yield when CO2 decearse molecules/J
	//double NO_net = NO_energy_yield_net*light_energy;//NO net increase/yr
	double NO_net = 2.6e9*30/14*0.066;// g yr-1,5e9 的单位是 g N yr-1
	double amount_net_NO = NO_net ;
	/*nitrogen fixation by methanogens*/
	//double n_N2;// nitrogen consumption of nitrogen fixation by methanogens, unit in mol
	//double n_ATP;// energy consumption of nitrogen fixation by methanogens, unit in mol
	//double n_H2;// hydrogen production of nitrogen fixation by methanogens, unit in mol




	int Hr = 24;
	double day_temp = Hr / 24;
	//定义CO2，H2和CO2从分压向溶解度的转换

	const double Hh2 = 0.00078*3, Hch4 = 0.0014*3;//define henry constant of h2 and ch4,(unit in mole/L*atm)
	float cH2, cCO2, cCH4;                  //define concentration of H2, CO2, CH4, (unit in mole/L)
	

	//通过浓度计算Q
	float Q;
	
	//通过Q计算detaG1
	//int detaG0 = -131;          //define the standard change of gibbs energy(unit in KJ/mole)
	double detaG0 = (-253 + 0.41*T);//KJ/mole Kral 1999
	const double R = 0.008314;  //define the gas constant(unit in KJ/mole/K)
	double detaG1;               //define the change of gibbs energy
	
	

	//计算甲烷菌质量增长
	float td=0;// = log(2) / specific_growth_rate;
	//float day_times = Hr / td;
	//float biomass_times = pow(2, day_times);
	double mglobal = 2.57e17;// 1.2e15;//define the NPP by P limitations
	double mb , mb_temp;
	double mburial = 0;
	//mb = m0 * biomass_times; //计算甲烷菌终产量
	
	double m0_temp = m0 ;

	/*质量回收*/
	//double m_de = 0.98*mb;//分解作用
	//计算K2
	
	//double rtns2K;//defiene how many times does reaction of methanogens biomass happen per unit time.(rnts/s)
	double Mmole = 23.28; // mole mass of methanogens biomass,(unit in g/mole)
	
	//double rest_m_N;//mount_ini_NO / 365 - 0.02* m_N + amount_net_NO / 365;//the rest N of the environment
	double nb, nH2, nCO2, VH2, VCO2, nCO2_de, VCO2_de, VCH4_loss;
	double m_sum = 0;
	double nCH4_m, nCH4 = 0, VCH4, nH2_m, nCO2_m, nCH4_loss = pCH4 / 10000 / 365, nH2fromCH4 = 3.5e12;//pCH4/10000/365;
	double mr_temp, mr=0;
	float specific_death_rate, specific_growth_rate_temp;
	double death_td, death_day_times, death_mass_times;
	int quzhengdetaday;
	double re_detaday = 0;// = quzhengdetaday % 365;
	float temperature;
	float growthrate;
	float Tdegree, Carbon, Hydrogen,Carbonindex;
	float Methane, Carbondioxide, Tedegree,Methaneindex;
	float T_h2_co2_min_plus, T_h2_co2_min_minus, T_h2_co2_max_plus, T_h2_co2_max_minus;
	float T_ch4_co2_min_plus, T_ch4_co2_min_minus, T_ch4_co2_max_plus, T_ch4_co2_max_minus;
	float T_h2_co2, T_ch4_co2;
	double minlenght, minlenght_h2, minlenght_temp, minlenght_temp_plus, minlenght_temp_minus, minlenghtch4, minlenghtch4_temp, minlenght1, minlenght1_temp = 1, minlenghtch4_temp_max=0.1, minlenghtch4_temp_min=0.1;
	double minlenght_temp_h2_min_plus, minlenght_temp_h2_min_minus, pH2_temp_min_length_plus, pH2_temp_min_length_minus, pH2_temp_max_length_minus, pH2_temp_max_length_plus, minlenght_temp_h2_max_plus, minlenght_temp_h2_max_minus;
	double pCH4_co2_temp_min_length_minus, pCH4_co2_temp_max_length_minus, pCH4_co2_temp_max_length_plus;
	double growthrate_temp_min_length = 0;
	int detaG2 = -30;
	double day_metha=0;
	double Vmole=0;
	double Va = 4e21/3;// unit in L,3倍大气压时大气体积
	double minlenght_ch4_co2, minlenghtch4_co2_temp_min_plus=0.1, pCH4_co2_temp_min_length_plus, minlenghtch4_co2_temp_min_minus=0.1, minlenghtch4_co2_temp_max_minus=0.1, minlenghtch4_co2_temp_max_plus=0.1;
	double distance_ph2_min_minus, distance_ph2_max_minus, abs_dis_ph2_min, abs_dis_ph2_max;
	double distance_ph2_min_plus, distance_ph2_max_plus;
	double T_h2_co2_min, T_h2_co2_max;
	double distance_pco2_min, distance_pco2_max, abs_dis_pco2;
	double distance_pch4_co2_min_minus, distance_pch4_co2_max_minus, abs_dis_pch4_co2_min_minus, abs_dis_pch4_co2_max_minus;
	double distance_pch4_co2_min_plus, distance_pch4_co2_max_plus, abs_dis_pch4_co2_min_plus, abs_dis_pch4_co2_max_plus;
	double abs_minlenghtch4;
	double pH2_temp_min_length = 0, pH2_temp_max_length, pCO2_temp_min_length = 0, pCO2_temp_max_length = 0;
	double pCH4_temp_min_length = 0, pCH4_temp_max_length = 0;
	double abs_minlenght;
	double T_ch4_co2_min, T_ch4_co2_max;
	double distance_pch4_min, distance_pch4_max, abs_dis_pch4;
	double abs_dis_ph2_min_minus, abs_dis_ph2_max_minus, abs_dis_ph2_min_plus, abs_dis_ph2_max_plus;
	double conH2 = 0, conCO2 = 0, conCH4 = 0, flux_h2 = 0, flux_ch4 = 0, m_h2, m_ch4 = 0;
	double vpx_H2 = 1.3e-2;//velocity of H2, unit in cm/s,*60*60*24  unit in cm/day
	double vpx_CH4 = 4.5e-3;//velocity of CH4,unit in cm/s
	double const Cindex = 6.02e20;//molecules cm-3 mol-1 L
	double nburial = mburial / Mmole;
	double ch4_OH_rate_const, ch4_O_rate_const;
	double conH2_temp, flux_h2_me=0, flux_ch4_me =0;
	double md, md_temp, mb_decompose, nCH4_dec, nCO2_dec, md_perday, n_nitrogen_dec=0;
	
	/*氢气部分*/
	ftot = pH2 + 2 * pCH4;//sum of partrial pressure
	H2_esc = 9e10 * ftot;//unit in hydrogen molecules/cm2/s
	/*last sentence， 9e10 should not be constant，it should be a function of ftot*/
	//H2_esc = 2.53e13 * ftot;//molecules/cm2/s
	//H2_esc = 4.9e11 * ftot;
	//H2_esc = 1e12 * ftot;
	H2_accum = H2_vol - H2_esc;//unit in hydrogen molecules/cm2/s
	//pH2_net = H2_accum / Na * 2 * S_at*t_year / m_at;//partial pressure by methanogens used 
	Vmole = (T / 273)*22.4 / 3;
	pH2_net = H2_accum / Na * S_at*t_year *Vmole/Va;
	//通过能量平衡计算K1
	cH2 = Hh2 * pH2 - H2_accum/Cindex/vpx_H2;                           //calculate the concentration of H2, CO2, CH4 by partial pressure and henry constant
	cCO2 = pCO2;
	cCH4 = pCH4*Hch4;// Va / Vmole / S_at / 10;//海洋表层的厚度为10m//20161122nCH4=pCH4*Va/Vmole
	//flux_ch4 = nCH4*Na / S_at / 60 / 60 / 24;
	//pCH4 = (flux_ch4 / Cindex / vpx_CH4 + cCH4) / Hch4;

	Q = cCH4 / cCO2 / pow(cH2, 4);
	detaG1 = detaG0 + R*T*log(Q);//calculate the change of gibbs energy by formula

	/*二氧化碳部分*/

	double Fsio2w = 3.9e11*pow(pCO2 / 0.00037, 0.3)*pow(2.7323, (T - 285) / 13.7); //mole/year
	double CO2_accumulate = Fridge + Farc - Fsio2w - Fhydro / 10;
	//double pCO2_net = CO2_accumulate * 44 / m_at;//methanogens cost--Fbio 
	double pCO2_net = CO2_accumulate * Vmole/Va;
	double m_N;
	//minlenghtch4_temp = 0.1*0.1 + 0.01*0.01;

	/*甲烷光解， Pavlov中 厌氧条件下甲烷的lifetime比现在长1000倍，现在的甲烷lifetime是10years，World Meteorological Organization. (1999)，所以早期地球甲烷的lifetime是10^4year,这个值在kasting 83年的文章中有提过*/
	if (detaG1 < 0)
	{
		
		double m_day;
		//float biomass_times = log10(mb / m0) / log10(2);//质量倍增了多少次
		td = log(2) / specific_growth_rate;//倍增一次的时间
		double times_day_biomass = 24 / td;//20161125
		mb = m0*pow(2, times_day_biomass);
		mburial = mb *0.02;//sensitivity test 20151231
		double nburial = mburial / Mmole;
		mb_temp = mb;
		m0_temp = m0;
		 m_N = mb * 14 / Mmole*0.24;//生物量中所包含的N;// = mb * 14 / Mmole*0.24;//生物量中所包含的N
		/*自然固氮部分*/
		//double rest_m_N = rest_m_N + amount_ini_NO - 0.02* m_N + amount_net_NO / 365;//the rest N of the environment
		//rest_m_N = rest_m_N + amount_ini_NO - 0.02* m_N - m_N + amount_net_NO / 365;//the rest N of the environment,0.002是埋藏率
		 
		m0 = mb;
		//detaday = detaday + 5;
		detaday++;
		year = detaday / 365.00;
		//nCH4_m = detaG1 / detaG0;
		//detaG2 = -1.36 * 30;

		if (m_N <= rest_m_N)
		{
			detaG2 = -1.36 * 30;// the energy for cell production(unit in KJ/mole) 

		}
		if (rest_m_N > 0 && m_N > rest_m_N)
		{
			double ratio_naturalnitrogenf_methanonitrogenf = rest_m_N / m_N;
			detaG2 = ratio_naturalnitrogenf_methanonitrogenf * (-1.36 * 30) + (1 - ratio_naturalnitrogenf_methanonitrogenf)*(-1.36 * 30 - 0.12 * 16 * 30);
		}
		if (rest_m_N <= 0)
		{
			detaG2 = -(1.36 * 30 + 0.12 * 16 * 30);
			rest_m_N = 0;
		}
		//times++;
		
		nCH4 = detaG2 / detaG1*mb / Mmole;
		nH2 = 4 * nCH4 + 2.09*mb / Mmole;
		nCO2 = nCH4 + mb / Mmole;
		Vmole = (T / 273)*22.4 / 3;
		ftot = pH2 + 2 * pCH4;
		//H2_esc = 2.53e13 * ftot;
		H2_esc = 9e10 * ftot;
		H2_accum = H2_vol - H2_esc;
		pH2_net = H2_accum / Na * S_at*t_year *Vmole / Va;//
		flux_ch4_me = nCH4*Na / S_at / 60 / 60 / 24;
		flux_h2_me = nH2*Na / S_at / 60 / 60 / 24;
		nCH4_loss = pCH4 / 10000 / 365;
		pCH4 = pCH4 + nCH4*Vmole/Va - nCH4_loss;//(cCH4 - flux_ch4_me / Cindex / vpx_CH4) / Hch4
		pH2 = pH2 - nH2*Vmole / Va + pH2_net / 365 + 0.5*nCH4_loss;
		pCO2 = pCO2 - nCO2*Vmole / Va + pCO2_net / 365 ;

		if (pCH4 < 1e-6)
		{
			pCH4 = 1e-6;
		}

		//conH2 = nH2 / S_at / 1;
		cH2 = Hh2 * pH2 - H2_accum / Cindex / vpx_H2;                           //calculate the concentration of H2, CO2, CH4 by partial pressure and henry constant
		cCO2 = pCO2;
		flux_ch4 = pCH4*Va / Vmole * Na / S_at / 60 / 60 / 24;
		cCH4 = fabs(flux_ch4 / vpx_CH4 / Cindex - Hch4*pCH4);
		
		rest_m_N = rest_m_N - m_N + amount_net_NO / 365 + n_nitrogen_dec;

		minlenght_temp_plus = 0.04;
		minlenght_temp_minus = 0.04;
		//minlenght_temp_h2_plus, minlenght_temp_h2_minus
		minlenght_temp_h2_min_plus = 0.11;
		minlenght_temp_h2_min_minus = 0.11;
		minlenght_temp_h2_max_plus = 0.11;
		minlenght_temp_h2_max_minus = 0.11;
		minlenghtch4_temp_max = 0.11;
		minlenghtch4_temp_min = 0.11;
		minlenghtch4_co2_temp_min_plus = 0.11;
		minlenghtch4_co2_temp_max_plus = 0.11;
		minlenghtch4_co2_temp_min_minus = 0.11;
		minlenghtch4_co2_temp_max_minus = 0.11;

		FILE *ww = fopen("pco2.txt", "r");//co2 data from Worthwords,20151222
		while (1)
		{
			//fscanf(ww, "%f %f %f", &Tdegree, &Carbon, &Hydrogen);
			//printf("%f %f %f\n",Tdegree,Carbon,Hydrogen);
			fscanf(ww, "%f", &Carbon);
			//printf("%f\n", Carbon);
			if (feof(ww))
			{
				break;
			}

			minlenght = pCO2 - Carbon;
			if (minlenght > 0)
			{
				if (minlenght <= minlenght_temp_plus)
				{
					minlenght_temp_plus = minlenght;
					//pH2_temp_min_length = Hydrogen;
					pCO2_temp_min_length = Carbon;
					//T_h2_co2 = Tdegree;
					//printf("%f\n", pCO2_temp_min_length);
				}
				else
				{
					continue;
				}
			}

			if (minlenght < 0)
			{
				abs_minlenght = fabs(minlenght);
				if (abs_minlenght <= minlenght_temp_minus)
				{
					minlenght_temp_minus = abs_minlenght;
					pCO2_temp_max_length = Carbon;

				}
				else
				{
					continue;
				}
			}

		}
		fclose(ww);
		FILE *ht = fopen("pco2_2.txt", "r");
		while (1)
		{
			fscanf(ht, "%f %f %f", &Carbonindex, &Hydrogen, &Tdegree);
			if (feof(ht))
			{
				break;
			}

			minlenght_h2 = pH2 - Hydrogen;
			if (minlenght_h2 > 0 && pCO2_temp_min_length == Carbonindex)
			{
				if (minlenght_h2 <= minlenght_temp_h2_min_plus)
				{
					minlenght_temp_h2_min_plus = minlenght_h2;
					pH2_temp_min_length_plus = Hydrogen;
					T_h2_co2_min_plus = Tdegree;
				}
				else
				{
					continue;
				}
			}
			if (minlenght_h2 < 0 && pCO2_temp_min_length == Carbonindex)
			{
				if (fabs(minlenght_h2) <= minlenght_temp_h2_min_minus)
				{
					minlenght_temp_h2_min_minus = fabs(minlenght_h2);
					pH2_temp_min_length_minus = Hydrogen;
					T_h2_co2_min_minus = Tdegree;
				}
				else
				{
					continue;
				}
			}
			if (minlenght_h2 < 0 && pCO2_temp_max_length == Carbonindex)
			{
				if (fabs(minlenght_h2) <= minlenght_temp_h2_max_minus)
				{
					minlenght_temp_h2_max_minus = fabs(minlenght_h2);
					pH2_temp_max_length_minus = Hydrogen;
					T_h2_co2_max_minus = Tdegree;
				}
				else
				{
					continue;
				}
			}
			if (minlenght_h2 > 0 && pCO2_temp_max_length == Carbonindex)
			{
				if (minlenght_h2 <= minlenght_temp_h2_max_plus)
				{
					minlenght_temp_h2_max_plus = minlenght_h2;
					pH2_temp_max_length_plus = Hydrogen;
					T_h2_co2_max_plus = Tdegree;
				}
				else
				{
					continue;
				}
			}

		}
		fclose(ht);
		distance_ph2_min_minus = pH2 - pH2_temp_min_length_minus;
		distance_ph2_min_plus = pH2 - pH2_temp_min_length_plus;
		distance_ph2_max_minus = pH2 - pH2_temp_max_length_minus;
		distance_ph2_max_plus = pH2 - pH2_temp_max_length_plus;
		abs_dis_ph2_min_minus = fabs(distance_ph2_min_minus);
		abs_dis_ph2_max_minus = fabs(distance_ph2_max_minus);
		abs_dis_ph2_min_plus = fabs(distance_ph2_min_plus);
		abs_dis_ph2_max_plus = fabs(distance_ph2_max_plus);
		T_h2_co2_min = abs_dis_ph2_min_plus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_minus + abs_dis_ph2_min_minus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_plus;
		T_h2_co2_max = abs_dis_ph2_max_plus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_minus + abs_dis_ph2_max_minus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_plus;
		distance_pco2_min = pCO2 - pCO2_temp_min_length;
		distance_pco2_max = pCO2 - pCO2_temp_max_length;
		abs_dis_pco2 = fabs(distance_pco2_max);
		T_h2_co2 = distance_pco2_min / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_max + abs_dis_pco2 / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_min;

		/*以下是考虑甲烷的温度*/
		FILE *me = fopen("pch4.txt", "r");
		while (1)
		{
			fscanf(me, "%f", &Methane);
			if (feof(me))
			{
				break;
			}
			minlenghtch4 = pCH4 - Methane;
			if (minlenghtch4 > 0)
			{
				if (minlenghtch4 <= minlenghtch4_temp_min)
				{
					minlenghtch4_temp_min = minlenghtch4;
					pCH4_temp_min_length = Methane;

				}

				else
				{
					continue;
				}
			}

			if (minlenghtch4 < 0)
			{
				abs_minlenghtch4 = fabs(minlenghtch4);
				if (abs_minlenghtch4 <= minlenghtch4_temp_max)
				{
					minlenghtch4_temp_max = abs_minlenghtch4;
					pCH4_temp_max_length = Methane;

				}

				else
				{
					continue;
				}
			}

		}
		fclose(me);

		FILE *hm = fopen("pch4_co2.txt", "r");//数据来自haqq-misra的文章，'a revised,hazy methane greenhouse for the archean earth',2008
		while (1)
		{
			fscanf(hm, "%f %f %f", &Methaneindex, &Carbondioxide, &Tedegree);
			if (feof(hm))
			{
				break;
			}

			minlenght_ch4_co2 = pCO2 - Carbondioxide;
			if (minlenght_ch4_co2 > 0 && pCH4_temp_min_length == Methaneindex)
			{
				if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_min_plus)
				{
					minlenghtch4_co2_temp_min_plus = minlenght_ch4_co2;
					pCH4_co2_temp_min_length_plus = Carbondioxide;
					T_ch4_co2_min_plus = Tedegree;

				}
				else
				{
					continue;
				}
			}
			if (minlenght_ch4_co2 < 0 && pCH4_temp_min_length == Methaneindex)
			{
				if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_min_minus)
				{
					minlenghtch4_co2_temp_min_minus = fabs(minlenght_ch4_co2);
					pCH4_co2_temp_min_length_minus = Carbondioxide;
					T_ch4_co2_min_minus = Tedegree;
				}
				else
				{
					continue;
				}
			}
			if (minlenght_ch4_co2 < 0 && pCH4_temp_max_length == Methaneindex)
			{
				if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_max_minus)
				{
					minlenghtch4_co2_temp_max_minus = fabs(minlenght_ch4_co2);
					pCH4_co2_temp_max_length_minus = Carbondioxide;
					T_ch4_co2_max_minus = Tedegree;
				}
				else
				{
					continue;
				}
			}
			if (minlenght_ch4_co2 > 0 && pCH4_temp_max_length == Methaneindex)
			{
				if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_max_plus)
				{
					minlenghtch4_co2_temp_max_plus = minlenght_ch4_co2;
					pCH4_co2_temp_max_length_plus = Carbondioxide;
					T_ch4_co2_max_plus = Tedegree;
				}
				else
				{
					continue;
				}
			}
		}
		fclose(hm);
		distance_pch4_co2_min_minus = pCO2 - pCH4_co2_temp_min_length_minus;
		distance_pch4_co2_min_plus = pCO2 - pCH4_co2_temp_min_length_plus;
		distance_pch4_co2_max_minus = pCO2 - pCH4_co2_temp_max_length_minus;
		distance_pch4_co2_max_plus = pCO2 - pCH4_co2_temp_max_length_plus;
		abs_dis_pch4_co2_min_minus = fabs(distance_pch4_co2_min_minus);
		abs_dis_pch4_co2_max_minus = fabs(distance_pch4_co2_max_minus);
		abs_dis_pch4_co2_min_plus = fabs(distance_pch4_co2_min_plus);
		abs_dis_pch4_co2_max_plus = fabs(distance_pch4_co2_max_plus);
		T_ch4_co2_min = abs_dis_pch4_co2_min_plus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_minus + abs_dis_pch4_co2_min_minus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_plus;
		T_ch4_co2_max = abs_dis_pch4_co2_max_plus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_minus + abs_dis_pch4_co2_max_minus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_plus;
		distance_pch4_min = pCH4 - pCH4_temp_min_length;
		distance_pch4_max = pCH4 - pCH4_temp_max_length;
		abs_dis_pch4 = fabs(distance_pch4_max);
		T_ch4_co2 = distance_pch4_min / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_max + abs_dis_pch4 / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_min;
		/*以下是温度的计算，考虑h2_co2和ch4_co2的共同作用*/

		T = (T_ch4_co2 + T_h2_co2) / 2;

		FILE *md_1 = fopen("a_b_inter.txt", "r");
		while (1)
		{
			fscanf(md_1, "%f %f\n", &temperature, &growthrate);
			double minlenght11 = temperature - T;
			minlenght1 = abs(minlenght11);

			if (minlenght1 < minlenght1_temp)
			{
				minlenght1_temp = minlenght1;
				growthrate_temp_min_length = temperature;
				specific_growth_rate = growthrate;
			}



			if (feof(md_1))
			{
				break;
			}
		}
		fclose(md_1);
		Q = cCH4 / cCO2 / pow(cH2, 4);
		detaG0 = (-253 + 0.41*T);
		detaG1 = detaG0 + R*T*log(Q);
		year = detaday / 365.00;
		//amount_ini_NO = 
		if (detaday == 1023)
		{
			printf("aaa");
		}
		FILE * yearcyclert;
		yearcyclert = fopen("365result20160220.txt", "a");
		//fprintf(yearcyclert, "%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
		//printf("%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
		fprintf(yearcyclert, " %d %f %f %f %f %f %f %f %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
		printf("%d %f %f %f %f %f %f %f  %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
		fclose(yearcyclert);

		
		if (H2_accum < flux_h2_me)
		{
			mr = mb;

			specific_growth_rate_temp = specific_growth_rate;
			//quzhengdetaday = int(detaday);
			//re_detaday = quzhengdetaday % 365;
			//double specific_growth_rate_temp, amount_ini_NO_temp;
			//specific_growth_rate = 1e-10;
			//specific_death_rate = specific_growth_rate_temp *0.5;
			//specific_death_rate = specific_growth_rate_temp * 0.25;
			//pecific_death_rate = 0.000408 ;
			//specific_death_rate = specific_growth_rate_temp * 0.1;//sensitivity analysis of death rate as follows
			//specific_death_rate = specific_growth_rate_temp * 0.05;
			//specific_death_rate = specific_growth_rate_temp * 0.01;
			//death_td = log(2) / specific_death_rate;
			//death_day_times = 24 / death_td;
			//death_mass_times = pow(0.5, death_day_times);
			while (H2_accum < flux_h2_me)
			{
				
				specific_growth_rate_temp = specific_growth_rate;
				specific_death_rate = specific_growth_rate_temp *0.5;
				death_td = log(2) / specific_death_rate;
				death_day_times = 24 / death_td;
				death_mass_times = pow(0.5, death_day_times);
				md = mr * death_mass_times;
				mr_temp = mr;
				if(md<1e-5)
				{
					md = 1e-5;
				}
				//mr_rest = mr *(1 - death_mass_times);
				mb = md;//mr_rest;
				mr = md;
				md_temp = md;

				md_perday = mr_temp - md_temp;
				nburial = md_perday*0.02/Mmole;
				mb_decompose = md_perday*0.98;
				n_nitrogen_dec = mb_decompose*0.24 * 14 / Mmole;
				//mb_decompose = 0.98*mb;
				//死亡部分分解产生的CO2 和 CH4 201611119
				nCH4_dec = mb_decompose/Mmole*0.5;
				nCO2_dec = mb_decompose/Mmole*0.5;
			
				//nH2 = 2.09*mb / Mmole;
				//VH2 = nH2*Vmole;
				//detaday = detaday+5;
				detaday++;
				year = detaday / 365.0000000;

				ftot = pH2 + 2 * pCH4;
				//H2_esc = 2.53e13 * ftot;
				H2_esc = 9e10 * ftot;
				H2_accum = H2_vol - H2_esc;
				nCH4 = detaG2 / detaG1*mb / Mmole+nCH4_dec;
				nH2 = 4 * nCH4 + 2.09*mb / Mmole;//H2_accum / Na * S_at * 60 * 60 * 24;//H2扩散速率限制甲烷菌能够利用的
				nCO2 = nCH4 + mb / Mmole;
				
				Vmole = (T / 273)*22.4 / 3;
				pH2_net = H2_accum / Na * S_at*t_year *Vmole / Va;
				Fsio2w = 3.9e11*pow(pCO2 / 0.00037, 0.3)*pow(2.7323, (T - 285) / 13.7);
				CO2_accumulate = Fridge + Farc - Fsio2w - Fhydro / 10;
				pCO2_net = CO2_accumulate * Vmole / Va;

				flux_h2_me = nH2*Na / S_at / 60 / 60 / 24;
				//flux_ch4_me = nCH4*Na / S_at / 60 / 60 / 24;
				pCH4 = pCH4 + nCH4*Vmole/Va - nCH4_loss; //(cCH4 - flux_ch4_me / Cindex / vpx_CH4) / Hch4
				pH2 = pH2 - nH2*Vmole / Va + pH2_net / 365  + 0.5*nCH4_loss;
				pCO2 = pCO2 - nCO2*Vmole / Va + pCO2_net / 365 + nCO2_dec*Vmole/Va - nburial *Vmole/Va;
				if (pH2 < 0)
				{
					pH2 = 0;
				}

				cH2 = Hh2 * pH2 - H2_accum / Cindex / vpx_H2;                           //calculate the concentration of H2, CO2, CH4 by partial pressure and henry constant
				cCO2 = pCO2;
				flux_ch4 = pCH4*Va / Vmole * Na / S_at / 60 / 60 / 24;
				cCH4 = fabs(flux_ch4 / vpx_CH4 / Cindex - Hch4*pCH4);

				m_N = mb * 14 / Mmole*0.24;
				NO_net = 2.6e9 * 30 / 14 * 0.066;// g yr-1,5e9 的单位是 g N yr-1
				amount_net_NO = NO_net ;
				if (rest_m_N < 0)
				{
					rest_m_N = 0;
				}
				//rest_m_N = rest_m_N + amount_net_NO + 0.98 * m_N; //自然固氮累加
				rest_m_N = rest_m_N + amount_net_NO/365  + n_nitrogen_dec;
				if (detaday == 1023)
				{
					printf("aaa");
				}
				FILE * yearcyclert;
				yearcyclert = fopen("365result20160220.txt", "a");
				//fprintf(yearcyclert, "%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				//printf("%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				fprintf(yearcyclert, " %d %f %f %f %f %f %f %f %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				printf("%d %f %f %f %f %f %f %f  %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				fclose(yearcyclert);
				//minlenght_temp = pow(pCO2, 2) + pow(pH2, 2);
				//minlenghtch4_temp = pow(pCO2, 2) + pow(pCH4, 2);
				minlenght_temp_plus = 0.04;
				minlenght_temp_minus = 0.04;
				//minlenght_temp_h2_plus, minlenght_temp_h2_minus
				minlenght_temp_h2_min_plus = 0.11;
				minlenght_temp_h2_min_minus = 0.11;
				minlenght_temp_h2_max_plus = 0.11;
				minlenght_temp_h2_max_minus = 0.11;
				minlenghtch4_temp_max = 0.11;
				minlenghtch4_temp_min = 0.11;
				minlenghtch4_co2_temp_min_plus = 0.11;
				minlenghtch4_co2_temp_max_plus = 0.11;
				minlenghtch4_co2_temp_min_minus = 0.11;
				minlenghtch4_co2_temp_max_minus = 0.11;

				FILE *ww = fopen("pco2.txt", "r");//co2 data from Worthwords,20151222
				while (1)
				{
					//fscanf(ww, "%f %f %f", &Tdegree, &Carbon, &Hydrogen);
					//printf("%f %f %f\n",Tdegree,Carbon,Hydrogen);
					fscanf(ww, "%f", &Carbon);
					//printf("%f\n", Carbon);
					if (feof(ww))
					{
						break;
					}

					minlenght = pCO2 - Carbon;
					if (minlenght > 0)
					{
						if (minlenght <= minlenght_temp_plus)
						{
							minlenght_temp_plus = minlenght;
							//pH2_temp_min_length = Hydrogen;
							pCO2_temp_min_length = Carbon;
							//T_h2_co2 = Tdegree;
							//printf("%f\n", pCO2_temp_min_length);
						}
						else
						{
							continue;
						}
					}

					if (minlenght < 0)
					{
						abs_minlenght = fabs(minlenght);
						if (abs_minlenght <= minlenght_temp_minus)
						{
							minlenght_temp_minus = abs_minlenght;
							pCO2_temp_max_length = Carbon;

						}
						else
						{
							continue;
						}
					}

				}
				fclose(ww);
				FILE *ht = fopen("pco2_2.txt", "r");
				while (1)
				{
					fscanf(ht, "%f %f %f", &Carbonindex, &Hydrogen, &Tdegree);
					if (feof(ht))
					{
						break;
					}

					minlenght_h2 = pH2 - Hydrogen;
					if (minlenght_h2 > 0 && pCO2_temp_min_length == Carbonindex)
					{
						if (minlenght_h2 <= minlenght_temp_h2_min_plus)
						{
							minlenght_temp_h2_min_plus = minlenght_h2;
							pH2_temp_min_length_plus = Hydrogen;
							T_h2_co2_min_plus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 < 0 && pCO2_temp_min_length == Carbonindex)
					{
						if (fabs(minlenght_h2) <= minlenght_temp_h2_min_minus)
						{
							minlenght_temp_h2_min_minus = fabs(minlenght_h2);
							pH2_temp_min_length_minus = Hydrogen;
							T_h2_co2_min_minus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 < 0 && pCO2_temp_max_length == Carbonindex)
					{
						if (fabs(minlenght_h2) <= minlenght_temp_h2_max_minus)
						{
							minlenght_temp_h2_max_minus = fabs(minlenght_h2);
							pH2_temp_max_length_minus = Hydrogen;
							T_h2_co2_max_minus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 > 0 && pCO2_temp_max_length == Carbonindex)
					{
						if (minlenght_h2 <= minlenght_temp_h2_max_plus)
						{
							minlenght_temp_h2_max_plus = minlenght_h2;
							pH2_temp_max_length_plus = Hydrogen;
							T_h2_co2_max_plus = Tdegree;
						}
						else
						{
							continue;
						}
					}

				}
				fclose(ht);
				distance_ph2_min_minus = pH2 - pH2_temp_min_length_minus;
				distance_ph2_min_plus = pH2 - pH2_temp_min_length_plus;
				distance_ph2_max_minus = pH2 - pH2_temp_max_length_minus;
				distance_ph2_max_plus = pH2 - pH2_temp_max_length_plus;
				abs_dis_ph2_min_minus = fabs(distance_ph2_min_minus);
				abs_dis_ph2_max_minus = fabs(distance_ph2_max_minus);
				abs_dis_ph2_min_plus = fabs(distance_ph2_min_plus);
				abs_dis_ph2_max_plus = fabs(distance_ph2_max_plus);
				T_h2_co2_min = abs_dis_ph2_min_plus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_minus + abs_dis_ph2_min_minus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_plus;
				T_h2_co2_max = abs_dis_ph2_max_plus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_minus + abs_dis_ph2_max_minus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_plus;
				distance_pco2_min = pCO2 - pCO2_temp_min_length;
				distance_pco2_max = pCO2 - pCO2_temp_max_length;
				abs_dis_pco2 = fabs(distance_pco2_max);
				T_h2_co2 = distance_pco2_min / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_max + abs_dis_pco2 / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_min;

				/*以下是考虑甲烷的温度*/
				FILE *me = fopen("pch4.txt", "r");
				while (1)
				{
					fscanf(me, "%f", &Methane);
					if (feof(me))
					{
						break;
					}
					minlenghtch4 = pCH4 - Methane;
					if (minlenghtch4 > 0)
					{
						if (minlenghtch4 <= minlenghtch4_temp_min)
						{
							minlenghtch4_temp_min = minlenghtch4;
							pCH4_temp_min_length = Methane;

						}

						else
						{
							continue;
						}
					}

					if (minlenghtch4 < 0)
					{
						abs_minlenghtch4 = fabs(minlenghtch4);
						if (abs_minlenghtch4 <= minlenghtch4_temp_max)
						{
							minlenghtch4_temp_max = abs_minlenghtch4;
							pCH4_temp_max_length = Methane;

						}

						else
						{
							continue;
						}
					}

				}
				fclose(me);

				FILE *hm = fopen("pch4_co2.txt", "r");//数据来自haqq-misra的文章，'a revised,hazy methane greenhouse for the archean earth',2008
				while (1)
				{
					fscanf(hm, "%f %f %f", &Methaneindex, &Carbondioxide, &Tedegree);
					if (feof(hm))
					{
						break;
					}

					minlenght_ch4_co2 = pCO2 - Carbondioxide;
					if (minlenght_ch4_co2 > 0 && pCH4_temp_min_length == Methaneindex)
					{
						if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_min_plus)
						{
							minlenghtch4_co2_temp_min_plus = minlenght_ch4_co2;
							pCH4_co2_temp_min_length_plus = Carbondioxide;
							T_ch4_co2_min_plus = Tedegree;

						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 < 0 && pCH4_temp_min_length == Methaneindex)
					{
						if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_min_minus)
						{
							minlenghtch4_co2_temp_min_minus = fabs(minlenght_ch4_co2);
							pCH4_co2_temp_min_length_minus = Carbondioxide;
							T_ch4_co2_min_minus = Tedegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 < 0 && pCH4_temp_max_length == Methaneindex)
					{
						if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_max_minus)
						{
							minlenghtch4_co2_temp_max_minus = fabs(minlenght_ch4_co2);
							pCH4_co2_temp_max_length_minus = Carbondioxide;
							T_ch4_co2_max_minus = Tedegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 > 0 && pCH4_temp_max_length == Methaneindex)
					{
						if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_max_plus)
						{
							minlenghtch4_co2_temp_max_plus = minlenght_ch4_co2;
							pCH4_co2_temp_max_length_plus = Carbondioxide;
							T_ch4_co2_max_plus = Tedegree;
						}
						else
						{
							continue;
						}
					}
				}
				fclose(hm);
				distance_pch4_co2_min_minus = pCO2 - pCH4_co2_temp_min_length_minus;
				distance_pch4_co2_min_plus = pCO2 - pCH4_co2_temp_min_length_plus;
				distance_pch4_co2_max_minus = pCO2 - pCH4_co2_temp_max_length_minus;
				distance_pch4_co2_max_plus = pCO2 - pCH4_co2_temp_max_length_plus;
				abs_dis_pch4_co2_min_minus = fabs(distance_pch4_co2_min_minus);
				abs_dis_pch4_co2_max_minus = fabs(distance_pch4_co2_max_minus);
				abs_dis_pch4_co2_min_plus = fabs(distance_pch4_co2_min_plus);
				abs_dis_pch4_co2_max_plus = fabs(distance_pch4_co2_max_plus);
				T_ch4_co2_min = abs_dis_pch4_co2_min_plus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_minus + abs_dis_pch4_co2_min_minus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_plus;
				T_ch4_co2_max = abs_dis_pch4_co2_max_plus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_minus + abs_dis_pch4_co2_max_minus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_plus;
				distance_pch4_min = pCH4 - pCH4_temp_min_length;
				distance_pch4_max = pCH4 - pCH4_temp_max_length;
				abs_dis_pch4 = fabs(distance_pch4_max);
				T_ch4_co2 = distance_pch4_min / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_max + abs_dis_pch4 / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_min;
				/*以下是温度的计算，考虑h2_co2和ch4_co2的共同作用*/

				T = (T_ch4_co2 + T_h2_co2) / 2;

				FILE *md_1 = fopen("a_b_inter.txt", "r");
				while (1)
				{
					fscanf(md_1, "%f %f\n", &temperature, &growthrate);
					double minlenght11 = temperature - T;
					minlenght1 = abs(minlenght11);

					if (minlenght1 < minlenght1_temp)
					{
						minlenght1_temp = minlenght1;
						growthrate_temp_min_length = temperature;
						specific_growth_rate = growthrate;
					}



					if (feof(md_1))
					{
						break;
					}
				}
				fclose(md_1);
				year = detaday / 365.00;

				if (mb < 1e-5)
				{
					mb = 1e-5;
					m0 = mb;
				}
				if (mb > 1e-5)
				{
					mb = mr;
					m0 = mb;
				}

				Q = cCH4 / cCO2 / pow(cH2, 4);
				detaG0 = (-253 + 0.41*T);
				detaG1 = detaG0 + R*T*log(Q);
				year = detaday / 365.00;
	
			}


		}
		if (T < 276.04 && T > 276 && mb >=1e-5)
		{
			mr = mb;
			//detaday = detaday++;
			//specific_growth_rate_temp = specific_growth_rate;
			//quzhengdetaday = int(detaday);
			//re_detaday = quzhengdetaday % 365;
			//specific_death_rate = 0.000408;
			//specific_death_rate = specific_growth_rate* 0.5;
			//death_td = log(2) / specific_death_rate;
			//death_day_times = 24 / death_td;
			//death_mass_times = pow(0.5, death_day_times);
			//mburial = mburial + mb;
			while (T < 276.04 && T >276 && mb >=1e-5)
			{
				
				specific_growth_rate_temp = specific_growth_rate;
				specific_death_rate = specific_growth_rate_temp *0.5;
				death_td = log(2) / specific_death_rate;
				death_day_times = 24 / death_td;
				death_mass_times = pow(0.5, death_day_times);
				md = mr * death_mass_times;
				mr_temp = mr;
				//mr_rest = mr *(1 - death_mass_times);
				mb = md;//mr_rest;
				mr = md;
				md_temp = md;

				md_perday = mr_temp - md_temp;
				nburial = md_perday / Mmole;
				n_nitrogen_dec = mb_decompose*0.24 * 14 / Mmole;
				//死亡部分分解产生的CO2 和 CH4 201611119


				detaday++;// = detaday + 5;

				if (detaday == 1023)
				{
				printf("aaa");
				}
				year = detaday / 365.0000000;
				//温度小于276.04K，甲烷菌不产生甲烷
				nCH4 = 0;//10 * detaG2 / detaG1*mb / Mmole;
				nH2 = 0;// 4 * nCH4 + 2.09*mb / Mmole;
				nCO2 = 0;// nCH4 + mb / Mmole;
				//也不发生分解作用
				nCH4_dec = 0;
				nCO2_dec = 0;
				//生物量被埋藏了20161119
				

				ftot = pH2 + 2 * pCH4;//sum of partrial pressure
				//H2_esc = 2.53e13 * ftot;//unit in hydrogen molecules/cm2/s
				/*last sentence， 9e10 should not be constant，it should be a function of ftot*/
				//H2_esc = 2.5e13 * ftot;//molecules/cm2/s
				H2_esc = 9e10 * ftot;
				//H2_esc = 1e12 * ftot;
				H2_accum = H2_vol - H2_esc;//unit in hydrogen molecules/cm2/s
				Vmole = (T / 273)*22.4 / 3;
				pH2_net = H2_accum / Na * S_at*t_year *Vmole / Va;//partial pressure by methanogens used 
				
				//Fsio2w = 1.62e19*pow(pCO2, 0.5)*pow(2.7323, (T - 290) / 13.7);	//silicate weathering on continents--
				Fsio2w = 3.9e11*pow(pCO2 / 0.00037, 0.3)*pow(2.7323, (T - 285) / 13.7);
				CO2_accumulate = Fridge + Farc - Fsio2w - Fhydro / 10;//mole/year
				if (CO2_accumulate < 0)
				{
					printf("%f\n", CO2_accumulate);
				}
				pCO2_net = CO2_accumulate * Vmole / Va;//methanogens cost--Fbio
				pCH4 = pCH4 - nCH4_loss ;
				if (pCH4 < 1e-6)
				{
					pCH4 = 1e-6;
				}
				pH2 = pH2 + pH2_net / 365  + 0.5*nCH4_loss;
				pCO2 = pCO2 + pCO2_net / 365 - nburial*Vmole/Va;//这个地方要改，因为pCO2_net是随温度变化的

				cH2 = Hh2*pH2 - flux_h2 / vpx_H2 / Cindex;
				cCO2 = pCO2;
				flux_ch4 = pCH4*Va/Vmole * Na / S_at / 60 / 60 / 24;
				cCH4 = fabs(flux_ch4 / vpx_CH4 / Cindex - Hch4*pCH4);
					m_N = mb * 14 / Mmole*0.24;
				NO_net = 2.6e9 * 30 / 14 * 0.066;// g yr-1,5e9 的单位是 g N yr-1
				amount_net_NO = NO_net ;
				rest_m_N = rest_m_N + amount_net_NO /365 + n_nitrogen_dec;
				minlenght_temp_plus = 0.04;
				minlenght_temp_minus = 0.04;
				//minlenght_temp_h2_plus, minlenght_temp_h2_minus
				minlenght_temp_h2_min_plus = 0.11;
				minlenght_temp_h2_min_minus = 0.11;
				minlenght_temp_h2_max_plus = 0.11;
				minlenght_temp_h2_max_minus = 0.11;
				minlenghtch4_temp_max = 0.11;
				minlenghtch4_temp_min = 0.11;
				minlenghtch4_co2_temp_min_plus = 0.11;
				minlenghtch4_co2_temp_max_plus = 0.11;
				minlenghtch4_co2_temp_min_minus = 0.11;
				minlenghtch4_co2_temp_max_minus = 0.11;

				FILE *ww = fopen("pco2.txt", "r");//co2 data from Worthwords,20151222
				while (1)
				{
					//fscanf(ww, "%f %f %f", &Tdegree, &Carbon, &Hydrogen); llllllllllldddddd
					//printf("%f %f %f\n",Tdegree,Carbon,Hydrogen);
					fscanf(ww, "%f", &Carbon);
					//printf("%f\n", Carbon);
					if (feof(ww))
					{
						break;
					}

					minlenght = pCO2 - Carbon;
					if (minlenght > 0)
					{
						if (minlenght <= minlenght_temp_plus)
						{
							minlenght_temp_plus = minlenght;
							//pH2_temp_min_length = Hydrogen;
							pCO2_temp_min_length = Carbon;
							//T_h2_co2 = Tdegree;
							//printf("%f\n", pCO2_temp_min_length);
						}
						else
						{
							continue;
						}
					}

					if (minlenght < 0)
					{
						abs_minlenght = fabs(minlenght);
						if (abs_minlenght <= minlenght_temp_minus)
						{
							minlenght_temp_minus = abs_minlenght;
							pCO2_temp_max_length = Carbon;

						}
						else
						{
							continue;
						}
					}

				}
				fclose(ww);
				FILE *ht = fopen("pco2_2.txt", "r");
				while (1)
				{
					fscanf(ht, "%f %f %f", &Carbonindex, &Hydrogen, &Tdegree);
					if (feof(ht))
					{
						break;
					}

					minlenght_h2 = pH2 - Hydrogen;
					if (minlenght_h2 > 0 && pCO2_temp_min_length == Carbonindex)
					{
						if (minlenght_h2 <= minlenght_temp_h2_min_plus)
						{
							minlenght_temp_h2_min_plus = minlenght_h2;
							pH2_temp_min_length_plus = Hydrogen;
							T_h2_co2_min_plus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 < 0 && pCO2_temp_min_length == Carbonindex)
					{
						if (fabs(minlenght_h2) <= minlenght_temp_h2_min_minus)
						{
							minlenght_temp_h2_min_minus = fabs(minlenght_h2);
							pH2_temp_min_length_minus = Hydrogen;
							T_h2_co2_min_minus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 < 0 && pCO2_temp_max_length == Carbonindex)
					{
						if (fabs(minlenght_h2) <= minlenght_temp_h2_max_minus)
						{
							minlenght_temp_h2_max_minus = fabs(minlenght_h2);
							pH2_temp_max_length_minus = Hydrogen;
							T_h2_co2_max_minus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 > 0 && pCO2_temp_max_length == Carbonindex)
					{
						if (minlenght_h2 <= minlenght_temp_h2_max_plus)
						{
							minlenght_temp_h2_max_plus = minlenght_h2;
							pH2_temp_max_length_plus = Hydrogen;
							T_h2_co2_max_plus = Tdegree;
						}
						else
						{
							continue;
						}
					}

				}
				fclose(ht);
				distance_ph2_min_minus = pH2 - pH2_temp_min_length_minus;
				distance_ph2_min_plus = pH2 - pH2_temp_min_length_plus;
				distance_ph2_max_minus = pH2 - pH2_temp_max_length_minus;
				distance_ph2_max_plus = pH2 - pH2_temp_max_length_plus;
				abs_dis_ph2_min_minus = fabs(distance_ph2_min_minus);
				abs_dis_ph2_max_minus = fabs(distance_ph2_max_minus);
				abs_dis_ph2_min_plus = fabs(distance_ph2_min_plus);
				abs_dis_ph2_max_plus = fabs(distance_ph2_max_plus);
				T_h2_co2_min = abs_dis_ph2_min_plus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_minus + abs_dis_ph2_min_minus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_plus;
				T_h2_co2_max = abs_dis_ph2_max_plus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_minus + abs_dis_ph2_max_minus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_plus;
				distance_pco2_min = pCO2 - pCO2_temp_min_length;
				distance_pco2_max = pCO2 - pCO2_temp_max_length;
				abs_dis_pco2 = fabs(distance_pco2_max);
				T_h2_co2 = distance_pco2_min / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_max + abs_dis_pco2 / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_min;

				/*以下是考虑甲烷的温度*/
				FILE *me = fopen("pch4.txt", "r");
				while (1)
				{
					fscanf(me, "%f", &Methane);
					if (feof(me))
					{
						break;
					}
					minlenghtch4 = pCH4 - Methane;
					if (minlenghtch4 > 0)
					{
						if (minlenghtch4 <= minlenghtch4_temp_min)
						{
							minlenghtch4_temp_min = minlenghtch4;
							pCH4_temp_min_length = Methane;

						}

						else
						{
							continue;
						}
					}

					if (minlenghtch4 < 0)
					{
						abs_minlenghtch4 = fabs(minlenghtch4);
						if (abs_minlenghtch4 <= minlenghtch4_temp_max)
						{
							minlenghtch4_temp_max = abs_minlenghtch4;
							pCH4_temp_max_length = Methane;

						}

						else
						{
							continue;
						}
					}

				}
				fclose(me);

				FILE *hm = fopen("pch4_co2.txt", "r");//数据来自haqq-misra的文章，'a revised,hazy methane greenhouse for the archean earth',2008
				while (1)
				{
					fscanf(hm, "%f %f %f", &Methaneindex, &Carbondioxide, &Tedegree);
					if (feof(hm))
					{
						break;
					}

					minlenght_ch4_co2 = pCO2 - Carbondioxide;
					if (minlenght_ch4_co2 > 0 && pCH4_temp_min_length == Methaneindex)
					{
						if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_min_plus)
						{
							minlenghtch4_co2_temp_min_plus = minlenght_ch4_co2;
							pCH4_co2_temp_min_length_plus = Carbondioxide;
							T_ch4_co2_min_plus = Tedegree;

						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 < 0 && pCH4_temp_min_length == Methaneindex)
					{
						if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_min_minus)
						{
							minlenghtch4_co2_temp_min_minus = fabs(minlenght_ch4_co2);
							pCH4_co2_temp_min_length_minus = Carbondioxide;
							T_ch4_co2_min_minus = Tedegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 < 0 && pCH4_temp_max_length == Methaneindex)
					{
						if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_max_minus)
						{
							minlenghtch4_co2_temp_max_minus = fabs(minlenght_ch4_co2);
							pCH4_co2_temp_max_length_minus = Carbondioxide;
							T_ch4_co2_max_minus = Tedegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 > 0 && pCH4_temp_max_length == Methaneindex)
					{
						if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_max_plus)
						{
							minlenghtch4_co2_temp_max_plus = minlenght_ch4_co2;
							pCH4_co2_temp_max_length_plus = Carbondioxide;
							T_ch4_co2_max_plus = Tedegree;
						}
						else
						{
							continue;
						}
					}
				} 
				fclose(hm);
				distance_pch4_co2_min_minus = pCO2 - pCH4_co2_temp_min_length_minus;
				distance_pch4_co2_min_plus = pCO2 - pCH4_co2_temp_min_length_plus;
				distance_pch4_co2_max_minus = pCO2 - pCH4_co2_temp_max_length_minus;
				distance_pch4_co2_max_plus = pCO2 - pCH4_co2_temp_max_length_plus;
				abs_dis_pch4_co2_min_minus = fabs(distance_pch4_co2_min_minus);
				abs_dis_pch4_co2_max_minus = fabs(distance_pch4_co2_max_minus);
				abs_dis_pch4_co2_min_plus = fabs(distance_pch4_co2_min_plus);
				abs_dis_pch4_co2_max_plus = fabs(distance_pch4_co2_max_plus);
				T_ch4_co2_min = abs_dis_pch4_co2_min_plus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_minus + abs_dis_pch4_co2_min_minus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_plus;
				T_ch4_co2_max = abs_dis_pch4_co2_max_plus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_minus + abs_dis_pch4_co2_max_minus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_plus;
				distance_pch4_min = pCH4 - pCH4_temp_min_length;
				distance_pch4_max = pCH4 - pCH4_temp_max_length;
				abs_dis_pch4 = fabs(distance_pch4_max);
				T_ch4_co2 = distance_pch4_min / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_max + abs_dis_pch4 / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_min;
				/*以下是温度的计算，考虑h2_co2和ch4_co2的共同作用*/
				T = (T_ch4_co2 + T_h2_co2) / 2;

				cH2 = Hh2*pH2 - flux_h2 / vpx_H2 / Cindex;
				cCO2 = pCO2;
				flux_ch4 = pCH4*Va / Vmole * Na / S_at / 60 / 60 / 24;
				cCH4 = fabs(flux_ch4 / vpx_CH4 / Cindex - Hch4*pCH4);

				Q = cCH4 / cCO2 / pow(cH2, 4);
				detaG0 = (-253 + 0.41*T);
				detaG1 = detaG0 + R*T*log(Q);
				year = detaday / 365.00;

				if (mb < 1e-5)
				{
					mb = 1e-5;
					m0 = mb;
				}
				if (mb > 1e-5)
				{
					mb = mr;
					m0 = mb;
				}
				FILE * yearcyclert;
				yearcyclert = fopen("365result20160220.txt", "a");
				//fprintf(yearcyclert, "%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				//printf("%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				fprintf(yearcyclert, " %d %f %f %f %f %f %f %f %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				printf("%d %f %f %f %f %f %f %f  %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				fclose(yearcyclert);
			}
		}
		
		if (T <= 276)
		{
			mr = mb;

			//specific_growth_rate_temp = specific_growth_rate;
			//quzhengdetaday = int(detaday);
			//re_detaday = quzhengdetaday % 365;
			//double specific_growth_rate_temp, amount_ini_NO_temp;
			//specific_growth_rate = 1e-10;
			//specific_death_rate = specific_growth_rate_temp *0.5;
			//specific_death_rate = specific_growth_rate_temp * 0.25;
			//specific_death_rate = 0.000408 ;
			//specific_death_rate = specific_growth_rate_temp * 0.1;//sensitivity analysis of death rate as follows
			//specific_death_rate = specific_growth_rate_temp * 0.05;
			//specific_death_rate = specific_growth_rate_temp * 0.01;
			//death_td = log(2) / specific_death_rate;
			//death_day_times = 24 / death_td * 365;
			//death_mass_times = pow(0.5, death_day_times);
			//mburial = mburial + mb;
			while (T <= 276)
			{

				specific_growth_rate_temp = specific_growth_rate;
				specific_death_rate = specific_growth_rate_temp *0.5;
				death_td = log(2) / specific_death_rate;
				death_day_times = 24 / death_td;
				death_mass_times = pow(0.5, death_day_times);
				md = mr * death_mass_times;
				mr_temp = mr;
				//mr_rest = mr *(1 - death_mass_times);
				mb = md;//mr_rest;
				mr = md;
				md_temp = md;

				md_perday = mr_temp - md_temp;
				nburial = md_perday/Mmole;
				

				detaday = detaday + 365;

				year = detaday / 365.0000000;
				
				nCH4 = 0;//10 * detaG2 / detaG1*mb / Mmole;
				nH2 = 0;// 4 * nCH4 + 2.09*mb / Mmole;
				nCO2 = 0;// nCH4 + mb / Mmole;
				//也不发生分解作用
				nCH4_dec = 0;
				nCO2_dec = 0;
				//生物量被埋藏了20161119
				//Fsio2w = 1.62e19*pow(pCO2, 0.5)*pow(2.7323, (T - 290) / 13.7);	//silicate weathering on continents--
				Fsio2w = 3.9e11*pow(pCO2 / 0.00037, 0.3)*pow(2.7323, (T - 285) / 13.7);
				CO2_accumulate = Fridge + Farc - Fsio2w - Fhydro / 10;
				if (CO2_accumulate < 0)
				{
					printf("%f\n", CO2_accumulate);
				}

				ftot = pH2 + 2 * pCH4;
				//H2_esc = 2.53e13 * ftot;
				H2_esc = 9e10 * ftot;
				//H2_esc = 1e12 * ftot;
				H2_accum = H2_vol - H2_esc;//unit in hydrogen molecules/cm2/s
				Vmole = (T / 273)*22.4 / 3;
				pH2_net = H2_accum / Na * S_at*t_year *Vmole / Va;
				Fsio2w = 3.9e11*pow(pCO2 / 0.00037, 0.3)*pow(2.7323, (T - 285) / 13.7);
				CO2_accumulate = Fridge + Farc - Fsio2w - Fhydro / 10;
				pCO2_net = CO2_accumulate * Vmole / Va;
				pCH4 = pCH4 - nCH4_loss * 365;
				if (pCH4 < 1e-6)
				{
					pCH4 = 1e-6;
				}
				pH2 = pH2 + pH2_net + 0.5*nCH4_loss * 365;
				pCO2 = pCO2 + pCO2_net;


				//flux_h2 = 0;
				
				flux_h2 = pH2*Va / Vmole*Na / S_at / 60 / 60 / 24;
				flux_h2_me = nH2*Na/S_at/60/60/24; //molecules cm-2 s-1
				flux_ch4_me = nCH4*Na / S_at / 60 / 60 / 24;
				//flux_ch4_me = 0;//molecules cm-2 s-1
				conH2 = fabs(Hh2*pH2 - flux_h2 / vpx_H2 / Cindex);
				cH2 = conH2;
				//cH2 = Hh2*pH2;
				cCO2 = pCO2;
				//conCH4 = fabs(flux_ch4_me / vpx_CH4 / Cindex - Hch4*pCH4);
				//cCH4 = conCH4;
				flux_ch4 = pCH4*Va / Vmole * Na / S_at / 60 / 60 / 24;
				cCH4 = fabs(flux_ch4 / vpx_CH4 / Cindex - Hch4*pCH4);

				m_N = mb * 14 / Mmole*0.24;
				NO_net = 2.6e9 * 30 / 14 * 0.066;// g yr-1,5e9 的单位是 g N yr-1
				amount_net_NO = NO_net;
				if (rest_m_N < 0)
				{
					rest_m_N = 0;
				}
				//rest_m_N = rest_m_N + amount_net_NO + 0.98 * m_N;//自然固氮累加
				rest_m_N = rest_m_N + amount_net_NO + n_nitrogen_dec;
				/*FILE * yearcyclert;
				yearcyclert = fopen("365result.txt", "a");
				//fprintf(yearcyclert, "%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				//printf("%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				fprintf(yearcyclert, " %d %f %d %f %f %f %f %f %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				printf("%d %f %d %f %f %f %f %f  %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				fclose(yearcyclert);*/
				//minlenght_temp = pow(pCO2, 2) + pow(pH2, 2);
				//minlenghtch4_temp = pow(pCO2, 2) + pow(pCH4, 2); 
				minlenght_temp_plus = 0.04;
				minlenght_temp_minus = 0.04;
				//minlenght_temp_h2_plus, minlenght_temp_h2_minus
				minlenght_temp_h2_min_plus = 0.11;
				minlenght_temp_h2_min_minus = 0.11;
				minlenght_temp_h2_max_plus = 0.11;
				minlenght_temp_h2_max_minus = 0.11;
				minlenghtch4_temp_max = 0.11;
				minlenghtch4_temp_min = 0.11;
				minlenghtch4_co2_temp_min_plus = 0.11;
				minlenghtch4_co2_temp_max_plus = 0.11;
				minlenghtch4_co2_temp_min_minus = 0.11;
				minlenghtch4_co2_temp_max_minus = 0.11;

				FILE *ww = fopen("pco2.txt", "r");//co2 data from Worthwords,20151222
				while (1)
				{
					//fscanf(ww, "%f %f %f", &Tdegree, &Carbon, &Hydrogen);
					//printf("%f %f %f\n",Tdegree,Carbon,Hydrogen);
					fscanf(ww, "%f", &Carbon);
					//printf("%f\n", Carbon);
					if (feof(ww))
					{
						break;
					}

					minlenght = pCO2 - Carbon;
					if (minlenght > 0)
					{
						if (minlenght <= minlenght_temp_plus)
						{
							minlenght_temp_plus = minlenght;
							//pH2_temp_min_length = Hydrogen;
							pCO2_temp_min_length = Carbon;
							//T_h2_co2 = Tdegree;
							//printf("%f\n", pCO2_temp_min_length);
						}
						else
						{
							continue;
						}
					}

					if (minlenght < 0)
					{
						abs_minlenght = fabs(minlenght);
						if (abs_minlenght <= minlenght_temp_minus)
						{
							minlenght_temp_minus = abs_minlenght;
							pCO2_temp_max_length = Carbon;

						}
						else
						{
							continue;
						}
					}

				}
				fclose(ww);
				FILE *ht = fopen("pco2_2.txt", "r");
				while (1)
				{
					fscanf(ht, "%f %f %f", &Carbonindex, &Hydrogen, &Tdegree);
					if (feof(ht))
					{
						break;
					}

					minlenght_h2 = pH2 - Hydrogen;
					if (minlenght_h2 > 0 && pCO2_temp_min_length == Carbonindex)
					{
						if (minlenght_h2 <= minlenght_temp_h2_min_plus)
						{
							minlenght_temp_h2_min_plus = minlenght_h2;
							pH2_temp_min_length_plus = Hydrogen;
							T_h2_co2_min_plus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 < 0 && pCO2_temp_min_length == Carbonindex)
					{
						if (fabs(minlenght_h2) <= minlenght_temp_h2_min_minus)
						{
							minlenght_temp_h2_min_minus = fabs(minlenght_h2);
							pH2_temp_min_length_minus = Hydrogen;
							T_h2_co2_min_minus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 < 0 && pCO2_temp_max_length == Carbonindex)
					{
						if (fabs(minlenght_h2) <= minlenght_temp_h2_max_minus)
						{
							minlenght_temp_h2_max_minus = fabs(minlenght_h2);
							pH2_temp_max_length_minus = Hydrogen;
							T_h2_co2_max_minus = Tdegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_h2 > 0 && pCO2_temp_max_length == Carbonindex)
					{
						if (minlenght_h2 <= minlenght_temp_h2_max_plus)
						{
							minlenght_temp_h2_max_plus = minlenght_h2;
							pH2_temp_max_length_plus = Hydrogen;
							T_h2_co2_max_plus = Tdegree;
						}
						else
						{
							continue;
						}
					}

				}
				fclose(ht);
				distance_ph2_min_minus = pH2 - pH2_temp_min_length_minus;
				distance_ph2_min_plus = pH2 - pH2_temp_min_length_plus;
				distance_ph2_max_minus = pH2 - pH2_temp_max_length_minus;
				distance_ph2_max_plus = pH2 - pH2_temp_max_length_plus;
				abs_dis_ph2_min_minus = fabs(distance_ph2_min_minus);
				abs_dis_ph2_max_minus = fabs(distance_ph2_max_minus);
				abs_dis_ph2_min_plus = fabs(distance_ph2_min_plus);
				abs_dis_ph2_max_plus = fabs(distance_ph2_max_plus);
				T_h2_co2_min = abs_dis_ph2_min_plus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_minus + abs_dis_ph2_min_minus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_plus;
				T_h2_co2_max = abs_dis_ph2_max_plus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_minus + abs_dis_ph2_max_minus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_plus;
				distance_pco2_min = pCO2 - pCO2_temp_min_length;
				distance_pco2_max = pCO2 - pCO2_temp_max_length;
				abs_dis_pco2 = fabs(distance_pco2_max);
				T_h2_co2 = distance_pco2_min / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_max + abs_dis_pco2 / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_min;

				/*以下是考虑甲烷的温度*/
				FILE *me = fopen("pch4.txt", "r");
				while (1)
				{
					fscanf(me, "%f", &Methane);
					if (feof(me))
					{
						break;
					}
					minlenghtch4 = pCH4 - Methane;
					if (minlenghtch4 > 0)
					{
						if (minlenghtch4 <= minlenghtch4_temp_min)
						{
							minlenghtch4_temp_min = minlenghtch4;
							pCH4_temp_min_length = Methane;

						}

						else
						{
							continue;
						}
					}

					if (minlenghtch4 < 0)
					{
						abs_minlenghtch4 = fabs(minlenghtch4);
						if (abs_minlenghtch4 <= minlenghtch4_temp_max)
						{
							minlenghtch4_temp_max = abs_minlenghtch4;
							pCH4_temp_max_length = Methane;

						}

						else
						{
							continue;
						}
					}

				}
				fclose(me);

				FILE *hm = fopen("pch4_co2.txt", "r");//数据来自haqq-misra的文章，'a revised,hazy methane greenhouse for the archean earth',2008
				while (1)
				{
					fscanf(hm, "%f %f %f", &Methaneindex, &Carbondioxide, &Tedegree);
					if (feof(hm))
					{
						break;
					}

					minlenght_ch4_co2 = pCO2 - Carbondioxide;
					if (minlenght_ch4_co2 > 0 && pCH4_temp_min_length == Methaneindex)
					{
						if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_min_plus)
						{
							minlenghtch4_co2_temp_min_plus = minlenght_ch4_co2;
							pCH4_co2_temp_min_length_plus = Carbondioxide;
							T_ch4_co2_min_plus = Tedegree;

						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 < 0 && pCH4_temp_min_length == Methaneindex)
					{
						if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_min_minus)
						{
							minlenghtch4_co2_temp_min_minus = fabs(minlenght_ch4_co2);
							pCH4_co2_temp_min_length_minus = Carbondioxide;
							T_ch4_co2_min_minus = Tedegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 < 0 && pCH4_temp_max_length == Methaneindex)
					{
						if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_max_minus)
						{
							minlenghtch4_co2_temp_max_minus = fabs(minlenght_ch4_co2);
							pCH4_co2_temp_max_length_minus = Carbondioxide;
							T_ch4_co2_max_minus = Tedegree;
						}
						else
						{
							continue;
						}
					}
					if (minlenght_ch4_co2 > 0 && pCH4_temp_max_length == Methaneindex)
					{
						if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_max_plus)
						{
							minlenghtch4_co2_temp_max_plus = minlenght_ch4_co2;
							pCH4_co2_temp_max_length_plus = Carbondioxide;
							T_ch4_co2_max_plus = Tedegree;
						}
						else
						{
							continue;
						}
					}
				}
				fclose(hm);
				distance_pch4_co2_min_minus = pCO2 - pCH4_co2_temp_min_length_minus;
				distance_pch4_co2_min_plus = pCO2 - pCH4_co2_temp_min_length_plus;
				distance_pch4_co2_max_minus = pCO2 - pCH4_co2_temp_max_length_minus;
				distance_pch4_co2_max_plus = pCO2 - pCH4_co2_temp_max_length_plus;
				abs_dis_pch4_co2_min_minus = fabs(distance_pch4_co2_min_minus);
				abs_dis_pch4_co2_max_minus = fabs(distance_pch4_co2_max_minus);
				abs_dis_pch4_co2_min_plus = fabs(distance_pch4_co2_min_plus);
				abs_dis_pch4_co2_max_plus = fabs(distance_pch4_co2_max_plus);
				T_ch4_co2_min = abs_dis_pch4_co2_min_plus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_minus + abs_dis_pch4_co2_min_minus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_plus;
				T_ch4_co2_max = abs_dis_pch4_co2_max_plus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_minus + abs_dis_pch4_co2_max_minus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_plus;
				distance_pch4_min = pCH4 - pCH4_temp_min_length;
				distance_pch4_max = pCH4 - pCH4_temp_max_length;
				abs_dis_pch4 = fabs(distance_pch4_max);
				T_ch4_co2 = distance_pch4_min / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_max + abs_dis_pch4 / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_min;
				/*以下是温度的计算，考虑h2_co2和ch4_co2的共同作用*/
				T = (T_ch4_co2 + T_h2_co2) / 2;

				Q = cCH4 / cCO2 / pow(cH2, 4);
				detaG0 = (-253 + 0.41*T);
				detaG1 = detaG0 + R*T*log(Q);
				year = detaday / 365.00;

				if (mb < 1e-5)
				{
					mb = 1e-5;
					m0 = mb;
				}
				if (mb > 1e-5)
				{
					mb = mr;
					m0 = mb;
				}
				FILE * yearcyclert;
				yearcyclert = fopen("365result20160220.txt", "a");
				//fprintf(yearcyclert, "%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				//printf("%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
				fprintf(yearcyclert, " %d %f %f %f %f %f %f %f %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				printf("%d %f %f %f %f %f %f %f  %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
				fclose(yearcyclert);
				if (runtimes == 5030929)
				{
					printf("aaa");
				}
			}


			//mb = m0;

			//year++;
			//detaday++;

		}

		

	

	}
		if (detaG1 >= 0)
		{
			//mb = m0;
			mr = mb;

			specific_growth_rate_temp = specific_growth_rate;
			quzhengdetaday = int(detaday);
			re_detaday = quzhengdetaday % 365;
			//double specific_growth_rate_temp, amount_ini_NO_temp;
			//specific_growth_rate = 1e-10;
			//specific_death_rate = specific_growth_rate_temp *0.5;
			//specific_death_rate = specific_growth_rate_temp * 0.25;
			specific_death_rate = 0.000408 ;
			//specific_death_rate = specific_growth_rate_temp * 0.1;//sensitivity analysis of death rate as follows
			//specific_death_rate = specific_growth_rate_temp * 0.05;
			//specific_death_rate = specific_growth_rate_temp * 0.01;
			death_td = log(2) / specific_death_rate;
			death_day_times = 24 / death_td ;
			death_mass_times = pow(0.5, death_day_times);
			while (detaG1 >= 0)
			{

				 
				specific_growth_rate_temp = specific_growth_rate;
				specific_death_rate = specific_growth_rate_temp *0.5;
				death_td = log(2) / specific_death_rate;
				death_day_times = 24 / death_td;
				death_mass_times = pow(0.5, death_day_times);
				md = mr * death_mass_times;
				mr_temp = mr;
				//mr_rest = mr *(1 - death_mass_times);
				mb = md;//mr_rest;
				mr = md;
				md_temp = md;

				md_perday = mr_temp - md_temp;
				nburial = md_perday*0.02/Mmole;
				mb_decompose = md_perday*0.98;
				//mb_decompose = 0.98*mb;
				//死亡部分分解产生的CO2 和 CH4 201611119
				nCH4_dec = mb_decompose/Mmole*0.5;
				nCO2_dec = mb_decompose/Mmole*0.5;
				n_nitrogen_dec = mb_decompose*0.24 * 14 / Mmole;
				detaday = detaday ;
				year = detaday / 365.0000000;
				//没有能量，甲烷菌不产生甲烷，不消耗二氧化碳和氢气
				nCH4 = 0;//10 * detaG2 / detaG1*mb / Mmole;
				nH2 = 0;//4 * nCH4 + 2.09*mb / Mmole;
				nCO2 = 0;//nCH4 + mb / Mmole;


				ftot = pH2 + 2 * pCH4;
				H2_esc = 9e10 * ftot;
				H2_accum = H2_vol - H2_esc;//unit in hydrogen molecules/cm2/s
				Vmole = (T / 273)*22.4 / 3;
				pH2_net = H2_accum / Na * S_at*t_year *Vmole / Va;
				Fsio2w = 3.9e11*pow(pCO2 / 0.00037, 0.3)*pow(2.7323, (T - 285) / 13.7);
				CO2_accumulate = Fridge + Farc - Fsio2w - Fhydro / 10;
				pCO2_net = CO2_accumulate * Vmole / Va;
				
				pCH4 = pCH4 + nCH4*Vmole / Va - nCH4_loss + nCH4_dec*Vmole/Va;
				pH2 = pH2 - nH2*Vmole / Va + pH2_net / 365  + 0.5*nCH4_loss;
				pCO2 = pCO2 - nCO2*Vmole / Va + pCO2_net / 365  + nCO2_dec*Vmole/Va-nburial*Vmole/Va;
				if (pH2 < 0)
				{
					pH2 = 0;
				}
				//pCH4 = pCH4 - nCH4_loss;
				if (pCH4 < 1e-6)
				{
					pCH4 = 1e-6;
				}
				cH2 = Hh2*pH2 - flux_h2 / vpx_H2 / Cindex;
				cCO2 = pCO2;
				flux_ch4 = pCH4*Va / Vmole * Na / S_at / 60 / 60 / 24;
				cCH4 = fabs(flux_ch4 / vpx_CH4 / Cindex - Hch4*pCH4);
			
				m_N = mb * 14 / Mmole*0.24;
				NO_net = 2.6e9 * 30 / 14 * 0.066;// g yr-1,5e9 的单位是 g N yr-1
				amount_net_NO = NO_net;
				if (rest_m_N < 0)
				{
					rest_m_N = 0;
				}
				//rest_m_N = rest_m_N + amount_net_NO + 0.98 * m_N;; //自然固氮累加
				rest_m_N = rest_m_N + amount_net_NO / 365+ n_nitrogen_dec;
				minlenght_temp_plus = 0.04;
				minlenght_temp_minus = 0.04;
				//minlenght_temp_h2_plus, minlenght_temp_h2_minus
				minlenght_temp_h2_min_plus = 0.11;
				minlenght_temp_h2_min_minus = 0.11;
				minlenght_temp_h2_max_plus = 0.11;
				minlenght_temp_h2_max_minus = 0.11;
				minlenghtch4_temp_max = 0.11;
				minlenghtch4_temp_min = 0.11;
				minlenghtch4_co2_temp_min_plus = 0.11;
				minlenghtch4_co2_temp_max_plus = 0.11;
				minlenghtch4_co2_temp_min_minus = 0.11;
				minlenghtch4_co2_temp_max_minus = 0.11;

			FILE *ww = fopen("pco2.txt", "r");//co2 data from Worthwords,20151222
			while (1)
			{
				//fscanf(ww, "%f %f %f", &Tdegree, &Carbon, &Hydrogen);
				//printf("%f %f %f\n",Tdegree,Carbon,Hydrogen);
				fscanf(ww, "%f", &Carbon);
				//printf("%f\n", Carbon);
				if (feof(ww))
				{
					break;
				}

				minlenght = pCO2 - Carbon;
				if (minlenght > 0)
				{
					if (minlenght <= minlenght_temp_plus)
					{
						minlenght_temp_plus = minlenght;
						//pH2_temp_min_length = Hydrogen;
						pCO2_temp_min_length = Carbon;
						//T_h2_co2 = Tdegree;
						//printf("%f\n", pCO2_temp_min_length);
					}
					else
					{
						continue;
					}
				}

				if (minlenght < 0)
				{
					abs_minlenght = fabs(minlenght);
					if (abs_minlenght <= minlenght_temp_minus)
					{
						minlenght_temp_minus = abs_minlenght;
						pCO2_temp_max_length = Carbon;

					}
					else
					{
						continue;
					}
				}

			}
			fclose(ww);
			FILE *ht = fopen("pco2_2.txt", "r");
			while (1)
			{
				fscanf(ht, "%f %f %f", &Carbonindex, &Hydrogen, &Tdegree);
				if (feof(ht))
				{
					break;
				}

				minlenght_h2 = pH2 - Hydrogen;
				if (minlenght_h2 > 0 && pCO2_temp_min_length == Carbonindex)
				{
					if (minlenght_h2 <= minlenght_temp_h2_min_plus)
					{
						minlenght_temp_h2_min_plus = minlenght_h2;
						pH2_temp_min_length_plus = Hydrogen;
						T_h2_co2_min_plus = Tdegree;
					}
					else
					{
						continue;
					}
				}
				if (minlenght_h2 < 0 && pCO2_temp_min_length == Carbonindex)
				{
					if (fabs(minlenght_h2) <= minlenght_temp_h2_min_minus)
					{
						minlenght_temp_h2_min_minus = fabs(minlenght_h2);
						pH2_temp_min_length_minus = Hydrogen;
						T_h2_co2_min_minus = Tdegree;
					}
					else
					{
						continue;
					}
				}
				if (minlenght_h2 < 0 && pCO2_temp_max_length == Carbonindex)
				{
					if (fabs(minlenght_h2) <= minlenght_temp_h2_max_minus)
					{
						minlenght_temp_h2_max_minus = fabs(minlenght_h2);
						pH2_temp_max_length_minus = Hydrogen;
						T_h2_co2_max_minus = Tdegree;
					}
					else
					{
						continue;
					}
				}
				if (minlenght_h2 > 0 && pCO2_temp_max_length == Carbonindex)
				{
					if (minlenght_h2 <= minlenght_temp_h2_max_plus)
					{
						minlenght_temp_h2_max_plus = minlenght_h2;
						pH2_temp_max_length_plus = Hydrogen;
						T_h2_co2_max_plus = Tdegree;
					}
					else
					{
						continue;
					}
				}

			}
			fclose(ht);
			distance_ph2_min_minus = pH2 - pH2_temp_min_length_minus;
			distance_ph2_min_plus = pH2 - pH2_temp_min_length_plus;
			distance_ph2_max_minus = pH2 - pH2_temp_max_length_minus;
			distance_ph2_max_plus = pH2 - pH2_temp_max_length_plus;
			abs_dis_ph2_min_minus = fabs(distance_ph2_min_minus);
			abs_dis_ph2_max_minus = fabs(distance_ph2_max_minus);
			abs_dis_ph2_min_plus = fabs(distance_ph2_min_plus);
			abs_dis_ph2_max_plus = fabs(distance_ph2_max_plus);
			T_h2_co2_min = abs_dis_ph2_min_plus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_minus + abs_dis_ph2_min_minus / (abs_dis_ph2_min_minus + abs_dis_ph2_min_plus)*T_h2_co2_min_plus;
			T_h2_co2_max = abs_dis_ph2_max_plus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_minus + abs_dis_ph2_max_minus / (abs_dis_ph2_max_minus + abs_dis_ph2_max_plus)*T_h2_co2_max_plus;
			distance_pco2_min = pCO2 - pCO2_temp_min_length;
			distance_pco2_max = pCO2 - pCO2_temp_max_length;
			abs_dis_pco2 = fabs(distance_pco2_max);
			T_h2_co2 = distance_pco2_min / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_max + abs_dis_pco2 / (distance_pco2_min + abs_dis_pco2)*T_h2_co2_min;

			/*以下是考虑甲烷的温度*/
			FILE *me = fopen("pch4.txt", "r");
			while (1)
			{
				fscanf(me, "%f", &Methane);
				if (feof(me))
				{
					break;
				}
				minlenghtch4 = pCH4 - Methane;
				if (minlenghtch4 > 0)
				{
					if (minlenghtch4 <= minlenghtch4_temp_min)
					{
						minlenghtch4_temp_min = minlenghtch4;
						pCH4_temp_min_length = Methane;

					}

					else
					{
						continue;
					}
				}

				if (minlenghtch4 < 0)
				{
					abs_minlenghtch4 = fabs(minlenghtch4);
					if (abs_minlenghtch4 <= minlenghtch4_temp_max)
					{
						minlenghtch4_temp_max = abs_minlenghtch4;
						pCH4_temp_max_length = Methane;

					}

					else
					{
						continue;
					}
				}

			}
			fclose(me);

			FILE *hm = fopen("pch4_co2.txt", "r");//数据来自haqq-misra的文章，'a revised,hazy methane greenhouse for the archean earth',2008
			while (1)
			{
				fscanf(hm, "%f %f %f", &Methaneindex, &Carbondioxide, &Tedegree);
				if (feof(hm))
				{
					break;
				}

				minlenght_ch4_co2 = pCO2 - Carbondioxide;
				if (minlenght_ch4_co2 > 0 && pCH4_temp_min_length == Methaneindex)
				{
					if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_min_plus)
					{
						minlenghtch4_co2_temp_min_plus = minlenght_ch4_co2;
						pCH4_co2_temp_min_length_plus = Carbondioxide;
						T_ch4_co2_min_plus = Tedegree;

					}
					else
					{
						continue;
					}
				}
				if (minlenght_ch4_co2 < 0 && pCH4_temp_min_length == Methaneindex)
				{
					if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_min_minus)
					{
						minlenghtch4_co2_temp_min_minus = fabs(minlenght_ch4_co2);
						pCH4_co2_temp_min_length_minus = Carbondioxide;
						T_ch4_co2_min_minus = Tedegree;
					}
					else
					{
						continue;
					}
				}
				if (minlenght_ch4_co2 < 0 && pCH4_temp_max_length == Methaneindex)
				{
					if (fabs(minlenght_ch4_co2) <= minlenghtch4_co2_temp_max_minus)
					{
						minlenghtch4_co2_temp_max_minus = fabs(minlenght_ch4_co2);
						pCH4_co2_temp_max_length_minus = Carbondioxide;
						T_ch4_co2_max_minus = Tedegree;
					}
					else
					{
						continue;
					}
				}
				if (minlenght_ch4_co2 > 0 && pCH4_temp_max_length == Methaneindex)
				{
					if (minlenght_ch4_co2 <= minlenghtch4_co2_temp_max_plus)
					{
						minlenghtch4_co2_temp_max_plus = minlenght_ch4_co2;
						pCH4_co2_temp_max_length_plus = Carbondioxide;
						T_ch4_co2_max_plus = Tedegree;
					}
					else
					{
						continue;
					}
				}
			}
			fclose(hm);
			distance_pch4_co2_min_minus = pCO2 - pCH4_co2_temp_min_length_minus;
			distance_pch4_co2_min_plus = pCO2 - pCH4_co2_temp_min_length_plus;
			distance_pch4_co2_max_minus = pCO2 - pCH4_co2_temp_max_length_minus;
			distance_pch4_co2_max_plus = pCO2 - pCH4_co2_temp_max_length_plus;
			abs_dis_pch4_co2_min_minus = fabs(distance_pch4_co2_min_minus);
			abs_dis_pch4_co2_max_minus = fabs(distance_pch4_co2_max_minus);
			abs_dis_pch4_co2_min_plus = fabs(distance_pch4_co2_min_plus);
			abs_dis_pch4_co2_max_plus = fabs(distance_pch4_co2_max_plus);
			T_ch4_co2_min = abs_dis_pch4_co2_min_plus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_minus + abs_dis_pch4_co2_min_minus / (abs_dis_pch4_co2_min_minus + abs_dis_pch4_co2_min_plus)*T_ch4_co2_min_plus;
			T_ch4_co2_max = abs_dis_pch4_co2_max_plus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_minus + abs_dis_pch4_co2_max_minus / (abs_dis_pch4_co2_max_minus + abs_dis_pch4_co2_max_plus)*T_ch4_co2_max_plus;
			distance_pch4_min = pCH4 - pCH4_temp_min_length;
			distance_pch4_max = pCH4 - pCH4_temp_max_length;
			abs_dis_pch4 = fabs(distance_pch4_max);
			T_ch4_co2 = distance_pch4_min / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_max + abs_dis_pch4 / (distance_pch4_min + abs_dis_pch4)*T_ch4_co2_min;
			/*以下是温度的计算，考虑h2_co2和ch4_co2的共同作用*/
			T = (T_ch4_co2 + T_h2_co2) / 2;
			FILE *md_1 = fopen("a_b_inter.txt", "r");
			while (1)
			{
				fscanf(md_1, "%f %f\n", &temperature, &growthrate);
				double minlenght11 = temperature - T;
				minlenght1 = abs(minlenght11);

				if (minlenght1 < minlenght1_temp)
				{
					minlenght1_temp = minlenght1;
					growthrate_temp_min_length = temperature;
					specific_growth_rate = growthrate;
				}



				if (feof(md_1))
				{
					break;
				}
			}
			fclose(md_1);
			year = detaday / 365.00;

			if (mb < 1e-5)
			{
				mb = 1e-5;
				m0 = mb;
			}
			if (mb > 1e-5)
			{
				mb = mr;
				m0 = mb;
			}

			Q = cCH4 / cCO2 / pow(cH2, 4);
			detaG0 = (-253 + 0.41*T);
			detaG1 = detaG0 + R*T*log(Q);
			year = detaday / 365.0000000;
			if (runtimes == 5030929)
			{
				printf("aaa");
			}
			FILE * yearcyclert;
			yearcyclert = fopen("365result20160220.txt", "a");
			//fprintf(yearcyclert, "%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
			//printf("%f %f %d %f %f %f %f %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, mb_temp, m0_temp, m0, runtimes);
			fprintf(yearcyclert, " %d %f %f %f %f %f %f %f %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
			printf("%d %f %f %f %f %f %f %f  %f %d %f %f %d\n", detaday, year, T, pH2, pCO2, pCH4, m0_temp, mb, detaG1, detaG2, m_N, rest_m_N, runtimes);
			fclose(yearcyclert);

		}




	}


}

void main()
{

	//double pH2 = 0.1;
	
	double pH2 = 0.099990902;
	double pCO2 = 0.03070311;
	//double pCO2 = 0.03162;
	double pCH4 = 0.000264811;
	//double pCH4 = 1e-6;
	//double c_ace = 1e-10;//mole/L
	float specific_growth_rate = 0.009354213;
	double m0 = 1.52662E+11;//1e-5;
	float T = 283.114319;//281.793854;
	int detaday = 12677603;//16938670;//the change of day
	double year = 34733.1589;// 46407.315068;//total years
	int detayear = 0;//the change of year
	int times=0;
	//double rest_m_N = 2.6e9 * 30 / 14 * 1.5e9*0.066*1e-2;
	double rest_m_N = -15337827164;

	//double detaG1 = -131;

	double &rpH2 = pH2;
	double &rpCO2 = pCO2;
	double &rpCH4 = pCH4;
	float &rspecific_growth_rate = specific_growth_rate;
	double &rm0 = m0;
	float  &rT = T;
	int  &rdetaday = detaday;
	double  &ryear = year;//total years
	int  &rdetayear = detayear;//the change of year
	double &rrest_m_N = rest_m_N;
	//double &rdetaG1 = detaG1;

	int runtimes = 0;
	//while (runtimes != -1)	

	while(T>260)
	{ 
		runtimes++;
		calculate_Temperature(rpH2, rpCO2, rpCH4, rT, rspecific_growth_rate, rm0, rdetaday, ryear, rdetayear,rrest_m_N, times, runtimes);

		if (runtimes == 5030929)
		{
			printf("aaa");
		}

		
	}
	//return 0;
}
