#include "gmp.h"
#include <stdio.h>
#include <assert.h>
#include <string>
#include <cmath>
#include <iostream>
using namespace std;
void modular_addition(mpz_t a,mpz_t b,mpz_t c,mpz_t r)
{
	mpz_t x,y;
	mpz_init(x),mpz_init(y);
	mpz_mod(x,a,c);mpz_mod(y,b,c);
	mpz_add(r,x,y);
	mpz_mod(r,r,c);
	mpz_clear(x);mpz_clear(y);
}
void modular_multiplication(mpz_t a,mpz_t b,mpz_t c,mpz_t r)
{
	mpz_t x,y;
	mpz_init(x),mpz_init(y);
	mpz_mod(x,a,c);mpz_mod(y,b,c);
	mpz_mul(r,x,y);
	mpz_mod(r,r,c);
	mpz_clear(x);mpz_clear(y);
}
void modular_exponentiation(mpz_t a,mpz_t b,mpz_t c,mpz_t r)
{
	mpz_t p,tem;mpz_init(p);mpz_init(tem);
	if(mpz_cmp_ui(c,1)==0)
		mpz_set_ui(r,0);
	else
	{
		mpz_set_ui(r,1);
		mpz_mod(a,a,c);
		while(mpz_cmp_ui(b,0)>0)
		{
			mpz_mod_ui(p,b,2);
			if(mpz_cmp_ui(p,1)==0)
			{
				mpz_mul(tem,r,a);
				mpz_mod(r,tem,c);
			}
			mpz_fdiv_q_ui(b,b,2);
			mpz_mul(tem,a,a);
			mpz_mod(a,tem,c);
		}
	}
	mpz_clear(p);mpz_clear(tem);
}
void gcd(mpz_t a,mpz_t b)
{
	mpz_t t;mpz_init(t);
	while(mpz_cmp_ui(b,0)!=0)
	{
		mpz_set(t,b);
		mpz_mod(b,a,b);
		mpz_set(a,t);
	}
	mpz_clear(t);
}
void inverse_modulo(mpz_t r,mpz_t a,mpz_t c)
{
	mpz_t m0,t,q,x0,temp;
	mpz_init(m0);mpz_init(t);mpz_init(q);mpz_init(x0);mpz_init(temp);
	mpz_set(m0,c);mpz_set_ui(r,1);
	if(mpz_cmp_ui(c,1)==0)
		mpz_set_ui(r,0);
	while(mpz_cmp_ui(a,1)>0)
	{
		mpz_fdiv_q(q,a,c);
		mpz_set(t,c);
		mpz_mod(c,a,c);
		mpz_set(a,t);mpz_set(t,x0);
		mpz_mul(temp,q,x0);
		mpz_sub(x0,r,temp);
		mpz_set(r,t);
	}
	if(mpz_cmp_ui(r,0)<0)
		mpz_add(r,r,m0);
	mpz_clear(m0);mpz_clear(t);mpz_clear(q);mpz_clear(x0);mpz_clear(temp);
}
void baby_step_giant_step()
{
	yes:
	cout<<"\nEnter base, power and modulus [as in (base^x)=power(mod modulus)], where x is to be calculated:\n";
	mpz_t a,x,b,m;
	mpz_init(a);mpz_init(b);mpz_init(x);mpz_init(m);
	mpz_inp_str(a,stdin,10);
 	mpz_inp_str(b,stdin,10);
 	mpz_inp_str(m,stdin,10);
	mpz_t N,inv,w,GCD;long long int i,j;
	mpz_init(N);mpz_init(inv);mpz_init(w);mpz_init(GCD);
	mpz_gcd(GCD,a,m);
	if(mpz_cmp_ui(GCD,1)!=0)
	{
		cout<<"Sorry, base and modulus are not relatively prime.\n";
		cout<<"\nDo you want to do discrete logarithm again?  Press (y/n): ";char again;cin>>again;
		if(again=='y')
			goto yes;
		else
			return;
	}
	mpz_gcd(GCD,b,m);
	if(mpz_cmp_ui(GCD,1)!=0)
	{
		cout<<"Sorry, power and modulus are not relatively prime.\n";
		cout<<"\nDo you want to do discrete logarithm again?  Press (y/n): ";char again;cin>>again;
		if(again=='y')
			goto yes;
		else
			return;
	}
	if(mpz_cmp_ui(b,1)==0)
	{
		cout<<"Excluding zero, ";
	}
	mpz_sqrt(N,m);
	mpz_add_ui(N,N,1);
	mpz_mod(b,b,m);
	mpz_t *pow = new mpz_t[mpz_get_ui(N)-1];
	mpz_init(pow[0]);
	mpz_set(pow[0],a);
	for (i = 1 ; i < (mpz_get_ui (N)-1) ; i++) 
	{
      	mpz_init (pow[i]) ;
      	mpz_mul (pow[i], pow[i-1], a);
      	mpz_mod (pow[i], pow[i],m);
      	if(mpz_cmp(pow[i],b)==0)
      	{
      		cout<<"Required discrete logarithm = "<<i+1<<endl;
      		cout<<"NOTE: This is the smallest possible answer.\n";
      		return;
      	}
    }
    int count =0;
    mpz_set(w,N);
    for(i=0;i<=(mpz_get_ui(w)-2);i++)
    {
    	mpz_powm(inv,a,N,m);
    	mpz_invert(inv,inv,m);
    	mpz_mul(inv,inv,b);
    	mpz_mod(inv,inv,m);

    	for(j=0;j<=(mpz_get_ui(w)-2);j++)
    	{
    		if(mpz_cmp(pow[j],inv) ==0)
    		{
    			count=1;
    			mpz_add_ui(N,N,j+1);
    			mpz_set(x,N);
    			break;
    		}
    	}
    	if(count==1)
    	{
    		break;
    	}
    	mpz_add(N,N,w);
    }
    if(count==0)
    	cout<<"Sorry,discrete logarithm does not exist.\n";
    else
    {
    	cout<<"Required discrete logarithm = ";mpz_out_str(stdout,10,x);cout<<endl;
    	cout<<"NOTE: This is the smallest possible answer.\n";
    }
	mpz_clear(N);mpz_clear(inv);mpz_clear(w);mpz_clear(GCD);mpz_clear(a);mpz_clear(x);mpz_clear(b);mpz_clear(m);
}
void RSA_crpto()
{
	mpz_t p,q,n,e,phi,d,mes;
	mpz_init(p);mpz_init(q);mpz_init(n);mpz_init(e);mpz_init(phi);mpz_init(d);mpz_init(mes);
    unsigned int i, seed;
    gmp_randstate_t state;
    cout<<"Enter any random number (in the range of int) to generate RSA keys.\n";
	cin>>seed;
	cout<<":RSA key generation:\n";
	gmp_randinit_default(state);
    gmp_randseed_ui(state,seed);
	mpz_urandomb(p,state,512);
    mpz_urandomb(q,state,512);
    mpz_urandomb(e,state,17);
    mpz_nextprime(p,p);
    mpz_nextprime(q,q);
    mpz_nextprime(e,e);
    mpz_mul(n,p,q);
    cout<<"The generated RSA public key pair (n,e) is (";mpz_out_str(stdout,10,n);
    cout<<",";mpz_out_str(stdout,10,e);cout<<"), where n is modulus and e is public key exponent.\n";
    mpz_sub_ui(p,p,1);mpz_sub_ui(q,q,1);
    mpz_mul(phi,p,q);
	mpz_invert(d,e,phi);
	cout<<"\nThe generated RSA private key pair is ";mpz_out_str(stdout,10,d);cout<<".\n";
    cout<<"\n:RSA Encryption:\n( (m^e)mod n ) where m: message.";
    cout<<"\nEnter any number m (the message) that you want to encrypt.\n";
    mpz_inp_str(mes,stdin,10);
    cout<<"Encrypted message (cipher text): ";
    mpz_powm(mes,mes,e,n);
    mpz_out_str(stdout,10,mes);cout<<endl;
    cout<<"\n:RSA Decryption:";
    cout<<"\nDo you want to decrypt the above cipher text OR enter and decrypt a new cipher text of your own?  Enter 'a' for above OR 'n' for new.\n";
    char option;cin>>option;
    if(option=='a')
    {
    	cout<<"\nRSA Decryption. ( (c^d)mod n ) where c: cipher text.";
    	mpz_powm(mes,mes,d,n);
    	cout<<"\nDecrypted message (i.e the original message m): ";
    	mpz_out_str(stdout,10,mes);cout<<endl;
	}
	else
	{
		cout<<"Enter your cipher text.\n";
		mpz_inp_str(mes,stdin,10);
		cout<<"\n:RSA Decryption: ( (c^d)mod n ) where c: cipher text.";
    	mpz_powm(mes,mes,d,n);
    	cout<<"\nDecrypted message (i.e the original message m): ";
    	mpz_out_str(stdout,10,mes);cout<<endl;
	}
	gmp_randclear(state);
    mpz_clear(p);mpz_clear(q);mpz_clear(n);mpz_clear(e);mpz_clear(phi);mpz_clear(d);mpz_clear(mes);
}
void field_operations()
{
	cout<<"\nChoose your field.\n\t1. Prime fields\n\t2. Binary fields (only addition and subtraction)\nEnter your choice i.e. 1 or 2: ";
    int ch,check;cin>>ch;
    switch(ch)
    {
      case 1://prime fields
      cout<<"\nEnter your prime number p as in Fp(which denotes prime field).\n";
      mpz_t p,r;mpz_init(p);mpz_init(r);mpz_inp_str(p,stdin,10);
      check=mpz_probab_prime_p(p,40);
      if(check==0)
      {
      	  cout<<"Sorry, the number which you have input is not prime.";
      	  break;
      }
      cout<<"Select the operation that you want to perform.\n1. Addition\n2. Subtraction\n3. Multiplication\n4. Inversion\nEnter the corresponding number(i.e. 1 or 2 or 3 or 4): ";
      int ch1;cin>>ch1;
      switch(ch1)
      {
          	case 1://addition
          	cout<<"Enter two elements(numbers) of the prime field\n";
          	mpz_t a,b;mpz_init(a);mpz_init(b);
          	mpz_inp_str(a,stdin,10);mpz_inp_str(b,stdin,10);mpz_sub_ui(p,p,1);
          	if(mpz_cmp(a,p)>0||mpz_cmp(b,p)>0)
          	{
              	cout<<"Sorry,the inputted numbers do not belong to the field.\n";
              	break;
          	}
          	mpz_add_ui(p,p,1);
          	modular_addition(a,b,p,r);
          	mpz_out_str(stdout,10,a);cout<<" + ";mpz_out_str(stdout,10,b);cout<<" = ";mpz_out_str(stdout,10,r);
          	mpz_clear(a);mpz_clear(b);
          	break;
          	case 2://subtraction
          	cout<<"Enter two elements(numbers) of the prime field\n";
          	mpz_init(a);mpz_init(b);mpz_init(r);
          	mpz_inp_str(a,stdin,10);mpz_inp_str(b,stdin,10);mpz_sub_ui(p,p,1);
          	if(mpz_cmp(a,p)>0||mpz_cmp(b,p)>0)
          	{
             	cout<<"Sorry,the inputted numbers do not belong to the field.\n";
              	break;
          	}
          	mpz_add_ui(p,p,1);
          	if(mpz_cmp(a,b)>=0)
              	mpz_sub(r,a,b);
          	else
          	{
            	mpz_sub(r,a,b);
              	mpz_add(r,r,p);
          	}
          	mpz_out_str(stdout,10,a);cout<<" - ";mpz_out_str(stdout,10,b);cout<<" = ";mpz_out_str(stdout,10,r);
          	mpz_clear(a);mpz_clear(b);
          	break;
          	case 3://multiplication
          	cout<<"Enter two elements(numbers) of the prime field\n";
          	mpz_init(a);mpz_init(b);mpz_init(r);
          	mpz_inp_str(a,stdin,10);mpz_inp_str(b,stdin,10);mpz_sub_ui(p,p,1);
          	if(mpz_cmp(a,p)>0||mpz_cmp(b,p)>0)
          	{
              	cout<<"Sorry,the inputted numbers do not belong to the field.\n";
              	break;
          	}
          	mpz_add_ui(p,p,1);
          	modular_multiplication(a,b,p,r);
          	mpz_out_str(stdout,10,a);cout<<" x ";mpz_out_str(stdout,10,b);cout<<" = ";mpz_out_str(stdout,10,r);
          	mpz_clear(a);mpz_clear(b);
          	break;
          	case 4://inversion
          	mpz_init(a);mpz_init(r);
          	cout<<"Enter a number of the prime field which you want to invert.\n";
          	mpz_inp_str(a,stdin,10);mpz_sub_ui(p,p,1);
          	if(mpz_cmp(a,p)>0)
          	{
              	cout<<"Sorry,the inputted number does not belong to the field.\n";
              	break;
          	}
          	mpz_add_ui(p,p,1);
          	mpz_invert(r,a,p);
          	mpz_out_str(stdout,10,a);cout<<"^(-1) = ";mpz_out_str(stdout,10,r);
          	mpz_clear(a);
          	break;
          	default:
          	cout<<"Unacceptable choice.";
      	}
      	mpz_clear(p);mpz_clear(r);
      	break;
      	case 2://binary fields
      	cout<<"\nEnter the field size m of the field (as in Fp(2^m)): ";
      	int m;cin>>m;
      	cout<<"Select the operation that you want to perform.\n1. Addition\n2. Subtraction.\nEnter the corresponding number(i.e. 1 or 2): ";
        cin>>ch1;
        switch(ch1)
        {
        	case 1://addition
        	{cout<<"\nEnter two "<<m<<" bit strings representing corresponding polynomials.\n";
      		string s1,s2;cin>>s1>>s2;
      		cout<<s1<<" + "<<s2<<" = ";
      		for(int i=0;i<m;i++)
      		{
      			cout<<((s1[i]-'0')+(s2[i]-'0'))%2;
      		}cout<<endl;}
      		break;
      		case 2://subtraction
      		{cout<<"\nEnter two "<<m<<" bit strings representing corresponding polynomials.\n";
      		string s1,s2;cin>>s1>>s2;
      		cout<<s1<<" - "<<s2<<" = ";
      		for(int i=0;i<m;i++)
      		{
      			cout<<abs((s1[i]-'0')-(s2[i]-'0'));
      		}cout<<endl;}
      		break;
      		default:
      		cout<<"Unacceptable choice.";
      		return;
        }
      	break;
      	default:
      	cout<<"Unacceptable choice.";
    }
}
//main program
int main()
{
	int choice;
	jump:
	cout<<"\n\t\t::MENU::\n\t1. Modular Arithmetic\n\t2. GCD and inverse modulo n\n\t3. Discrete logarithm\n\t4. RSA Crptosystem\n\t5. Field operations in prime and binary fields.";
	cout<<"\nPlease enter your choice (i.e. 1 or 2 or 3 or 4 or 5): ";
	cin>>choice;
	switch(choice)
	{
		case 1://modular arithmetic
		cout<<"\nWhat do you want to do ?\n1. Modular Addition\n2. Modular Multiplication\n3. Modular Exponentiation";
		cout<<"\n\nEnter your choice (i.e. 1 or 2 or 3): ";
		int n1;cin>>n1;
		switch(n1)
		{
			case 1://modular addition
			mpz_t a,b,c,r;//r = result
  			mpz_init(a);mpz_init(b);mpz_init(c);mpz_init(r);
  			cout<<"\nEnter a,b and c [as in (a+b)mod c].\n";
  			mpz_inp_str(a,stdin,10);
 			mpz_inp_str(b,stdin,10);
 			mpz_inp_str(c,stdin,10);
 			modular_addition(a,b,c,r);
 			cout<<"(";mpz_out_str(stdout,10,a);cout<<" + ";mpz_out_str(stdout,10,b);cout<<") mod ";mpz_out_str(stdout,10,c);cout<<" = ";
 			mpz_out_str(stdout,10,r);
 			cout<<endl;
 			mpz_clear(a);mpz_clear(b);mpz_clear(c);mpz_clear(r);
 			cout<<"\nDo you want to do anything else? Enter (y/n): ";
 			char yn;cin>>yn;
 			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
			break;
			case 2://modular multiplication
			mpz_init(a);mpz_init(b);mpz_init(c);mpz_init(r);
  			cout<<"\nEnter a,b and c [as in (axb)mod c].\n";
  			mpz_inp_str(a,stdin,10);
 			mpz_inp_str(b,stdin,10);
 			mpz_inp_str(c,stdin,10);
 			modular_multiplication(a,b,c,r);
 			cout<<"(";mpz_out_str(stdout,10,a);cout<<" x ";mpz_out_str(stdout,10,b);cout<<") mod ";mpz_out_str(stdout,10,c);cout<<" = ";
 			mpz_out_str(stdout,10,r);
 			cout<<endl;
 			mpz_clear(a);mpz_clear(b);mpz_clear(c);mpz_clear(r);
 			cout<<"\nDo you want to do anything else? Enter (y/n): ";
 			cin>>yn;
 			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
  			break;
			case 3://modular exponentiation
			mpz_init(a);mpz_init(b);mpz_init(c);mpz_init(r);
  			cout<<"\nEnter base(a),exponent(b) and modulus(c) [as in (a^b)mod c].\n";
  			mpz_inp_str(a,stdin,10);
 			mpz_inp_str(b,stdin,10);
 			mpz_inp_str(c,stdin,10);
 			cout<<"(";mpz_out_str(stdout,10,a);cout<<" ^ ";mpz_out_str(stdout,10,b);cout<<") mod ";mpz_out_str(stdout,10,c);cout<<" = ";
 			modular_exponentiation(a,b,c,r);
 			mpz_out_str(stdout,10,r);
 			cout<<endl;
 			mpz_clear(a);mpz_clear(b);mpz_clear(c);mpz_clear(r);
 			cout<<"\nDo you want to do anything else? Enter (y/n): ";
 			cin>>yn;
 			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
			break;
			default:
			cout<<"Unacceptable choice\n";
			cout<<"Want to do again ? Press(y/n): ";
			cin>>yn;
			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
		}
		break;
		case 2://gcd and inverse modulo n
		cout<<"\nWhat do you want to do ?\n1. Find GCD of two numbers.\n2. Find inverse of a number (modulo n).";
		cout<<"\n\nEnter your choice (i.e. 1 or 2): ";
		int n2;cin>>n2;
		switch(n2)
		{
			case 1://gcd of two numbers
			mpz_t a,b;
			mpz_init(a);mpz_init(b);
 			cout<<"\nEnter the two numbers of whom you want to find GCD.\n";
  			mpz_inp_str(a,stdin,10);
 			mpz_inp_str(b,stdin,10);
			cout<<"GCD of ";mpz_out_str(stdout,10,a);cout<<" and ";mpz_out_str(stdout,10,b);cout<<" = ";
			gcd(a,b);
 			mpz_out_str(stdout,10,a);
 			cout<<endl;
 			mpz_clear(a);mpz_clear(b);
 			cout<<"\nDo you want to do anything else? Enter (y/n): ";
 			char yn;cin>>yn;
 			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
 			break;
 			case 2://inverse modulo n
 			mpz_t r,c,GCD;
 			mpz_init(a);mpz_init(c);mpz_init(r);mpz_init(GCD);
  			cout<<"\nEnter the number 'a' and the modulus 'n'.\n";
  			mpz_inp_str(a,stdin,10);
 			mpz_inp_str(c,stdin,10);
 			mpz_gcd(GCD,a,c);
 			if(mpz_cmp_ui(GCD,1) != 0)
 			{
 				cout<<"Sorry, inverse of ";mpz_out_str(stdout,10,a);cout<<" (modulo ";mpz_out_str(stdout,10,c);cout<<") does not exist since GCD OF ";
 				mpz_out_str(stdout,10,a);cout<<" and ";mpz_out_str(stdout,10,c);cout<<" IS NOT EQUAL TO 1.\n";
 				goto label;
 			}
 			cout<<"Inverse of ";mpz_out_str(stdout,10,a);cout<<" (modulo ";mpz_out_str(stdout,10,c);cout<<") = ";
 			inverse_modulo(r,a,c);
 			mpz_out_str(stdout,10,r);
 			cout<<endl;
 			label:
 			mpz_clear(a);mpz_clear(r);mpz_clear(c);mpz_clear(GCD);
 			cout<<"\nDo you want to do anything else? Enter (y/n): ";
 			cin>>yn;
 			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
			break;
 			default:
 			cout<<"\nUnacceptable Choice\n";
 			cout<<"Want to do again ? Press(y/n): ";
			cin>>yn;
			if(yn=='y')
 				goto jump;
 			else cout<<"THANK YOU !!\n\n";
 		}
		break;
		case 3://discrete logarithm
		cout<<"\nDiscrete logarithm.";
		baby_step_giant_step();
 		cout<<"\nDo you want to do anything else: Enter (y/n): ";
 		char yn;cin>>yn;
 		if(yn=='y')
 			goto jump;
 		else cout<<"THANK YOU !!\n\n";
		break;
		case 4://RSA Cryptosystem
		cout<<"\nRSA Cryptosystem.\n";
		RSA_crpto();
		cout<<"\nDo you want to do anything else? Enter (y/n): ";
 		cin>>yn;
 		if(yn=='y')
 			goto jump;
 		else cout<<"THANK YOU !!\n\n";
		break;
		case 5:
		field_operations();
		cout<<"\nDo you want to do anything else? Enter (y/n): ";
 		cin>>yn;
 		if(yn=='y')
 			goto jump;
 		else cout<<"THANK YOU !!\n\n";
		break;
		default:
		cout<<"\nUnacceptable Choice\n";
		cout<<"Want to do again ? Press(y/n): ";
		cin>>yn;
		if(yn=='y')
 			goto jump;
 		else cout<<"THANK YOU !!\n\n";
	}
	return 0;
}