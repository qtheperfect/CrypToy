#include <iostream>
#include <iomanip>
#include <ctime>
#include <cassert>
#include <functional>

// Requiring: libgmp-dev, libgmp10, libgmpxx4ldbl
#include <gmpxx.h> // g++ $filename.cpp -lgmp -lgmpxx

using namespace std;
typedef mpz_class bint;

//#define DEBUG

#define elapseTest(Exp) {                                       \
    cout<<"\n\n########## Begin speed test ########### \n\n";   \
    clock_t t1 = clock();                                       \
    Exp;                                                        \
    clock_t t2 = clock();                                       \
    double elapse = double(t2 - t1) / CLOCKS_PER_SEC ;          \
    cout<<"\n> Elapse (s): "<< setprecision(22)<<elapse <<endl; \
    cout<<"##########  End  speed test  ########## \n\n";       \
  }                                                             

namespace speeds {
  double funcElapse(function<void()> func) {
    double elapse = -1;
    cout<<"\n\n----------- Begin speed test -------------\n\n";
    clock_t t1 = clock();
    func();
    clock_t t2 = clock();
    elapse = double(t2 - t1) / CLOCKS_PER_SEC ;
    cout<<"\n> Elapse (s): "<< elapse <<endl;
    cout<<"\n\n-----------  End  speed test -------------\n\n";
    return elapse;
  }
}

namespace modcalcs_iter {

  bint modpow(bint base, bint exp, bint modulus) {
    bint result = 1;
    while (exp > 0) {
      if (exp % 2 == 0) {
        exp = exp / 2;
        base = base * base % modulus;
      }
      else {
        result = result * base % modulus;
        exp = (exp - 1) / 2;
        base = base * base % modulus;
      }
    }
    return result;
  }

  vector<bint> primeTesters = {2, 3, 11, 23, 47, 101};
  bool isPrime(bint n) {
    for (bint a : primeTesters) {
      if (modpow(a, n - 1, n) != 1)
        return false;
    }
    return true;
  }

  bint goodPrime12(bint n) {
    while (n > 11) {
      if ( isPrime(n) && isPrime((n - 1) / 2) )
        return n;
      else
        n -= 12;
    }
    return 11;
  }

  bint goodPrime(bint n) {
    return goodPrime12(n - (n % 12) + 11);
  }
  
}

namespace modcalcs {
  
  bint& modpow_ref(bint &base, bint &exp, bint &modulus, bint &result) {
    if (exp <= 0)
      return result;
    else if (exp % 2 == 0) {
      base = base * base % modulus;
      exp = exp / 2;
      return modpow_ref(base, exp, modulus, result);
    }
    else {
      result = result * base % modulus;
      exp = (exp - 1) / 2;
      base = base * base % modulus; 
      return modpow_ref(base, exp, modulus, result);
    }
  }
 
  bint modpow_const(const bint base, const bint exp, const bint modulus, const bint result) {
    //cout<<result<<" * "<<base<<" ^ "<<exp<<endl;
    if (exp <= 0)
      return result;
    else if (exp % 2 == 0) {
      return modpow_const(base * base % modulus, exp / 2, modulus, result);
    }
    else {
      //return modpow_const(base * base % modulus, (exp - 1) / 2, modulus, result * base % modulus);
      return modpow_const(base, exp - 1, modulus, result * base % modulus);
    }
  }
  
  
  bint modpow(bint base, bint exp, bint modulus) {
    bint result = bint(1);
    return modpow_ref(base, exp, modulus, result);
    //return modpow_const(base, exp, modulus, result);
  }

  bint& modpow(bint base, bint exp, bint modulus, bint &result) {
    return modpow_ref(base, exp, modulus, result);
  }

  vector<bint> primeTesters = {2, 3, 11, 23, 47, 101};
  bool isPrime(bint n) {
    for (bint a : primeTesters) {
      if (modpow(a, n - 1, n) != 1)
        return false;
    }
    return true;
  }

  bint& goodPrime12(bint& n) {
    if (isPrime(n) && isPrime((n - 1) / 2))
      return n;
    else if (n < 12) {
      n = 11;
      return n;
    }
    else {
      n = n - 12;
      return goodPrime12(n);
    }
  }

  bint goodPrime(bint n) {
    n -= (n % 12) - 11; 
    return goodPrime12(n);
  }

  class Bezouts {
  public:
    bint a, b;
    bint pa, pb, gcd;

    void getBezoutsRecur(bint& aa, bint& bb, bint& pa, bint& pb) {
      // find pa & pb so that pa * a + pb * b = gcd
      if (bb == 0) {
        pa = 1;
        this -> gcd = aa;

        // Assuming that &pb should be either & this->pa or & this->pb,
        if (& pb == & this->pb) 
          pb = 0;
        else
          pb = 1;
        // ... so that this->pa would always be positive,
        //     which is useful when seeking a modulo inverse of a.
        return ;
      }
      else {
        bint k = aa / bb;
        aa = aa % bb;
        getBezoutsRecur(bb, aa, pb, pa);
        pb = pb - k * pa;
        return ;
      }
    }
    void getBezoutsConst(const bint aa, const bint bb) {
      return getBezoutsConst(aa, bb, &this->pa, &this->pb);
    }
    void getBezoutsConst(const bint aa, const bint bb, bint* ap, bint* bp) {
      // find pa & pb so that pa * a + pb * b = gcd
      if (bb == 0) {
        *ap = 1;
        this -> gcd = aa;

        // Assuming that bp should be either & this->pa or & this->pb,
        if (bp == & this->pb) 
          *bp = 0;
        else
          *bp = 1;
        // ... so that this->pa would always be positive,
        //     which is useful when seeking a modulo inverse of a.
        return ;
      }
      else {
        bint k = aa / bb;
        getBezoutsConst(bb, aa % bb, bp, ap);
        *bp = *bp - k * *ap;
        return ;
      }
    }
    void getBezouts(bint aa, bint bb) {
      this -> a = aa;
      this -> b = bb;
      getBezoutsRecur(aa, bb, this -> pa, this -> pb);
      return;
    }
    void getBezouts(string aas, string bbs) {
      return getBezouts(bint(aas), bint(bbs));
    }
    void showFormula(){
      cout<<"(" << (pa) << ") * " << (a) << " + (" << (pb) << ") * " << b << " = " << ((pa * a + pb * b)) << " .... [gcd]-> " << (gcd)<<endl;
    }
    void showFormula(bint a, bint b, bint pa, bint pb){
      cout<<"(" << (pa) << ") * " << (a) << " + (" << (pb) << ") * " << b << " = " << ((pa * a + pb * b)) << " >-[gcd]-> " << (gcd)<<endl;
    }
    static void test() {
      cout<<" Beginning Bezout Transformation:\n\n";
      Bezouts bz = Bezouts();
      bz.getBezouts("243328423987719734823992491273832492847904", "483000001238947129347911200039289743200"); 
      bz.showFormula();
      bint modInverse(bint , bint);
      cout<<modInverse(bint("2384238428423783471"), bint("89347982379487972912412343"))<<endl;
      cout<<endl;
    }
  };

  bint modInverse(bint n, bint m) {
    Bezouts bz = Bezouts();
    bz.getBezouts(n, m);
    return bz.pa;
  }
}

namespace toyDSA {
  using namespace modcalcs;

  class AbstractDSA;
  class DSAPub;
  class DSA;
  class DSASign;

  class AbstractDSA {
  public:
    bint pubKey, q, p, g;
    AbstractDSA(){
      pubKey = q = p = g = 1;
    }
    AbstractDSA(bint pubKey, bint g, bint q, bint p) {
      this -> pubKey = pubKey ;
      this -> q = q;
      this -> p = p;
      this -> g = g;
    }
    bool check(const DSASign* sObj); 
  };

  class DSAPub : public AbstractDSA {
  public:
    DSAPub(bint pubKey, bint g, bint q, bint p):
      AbstractDSA(pubKey, g, q, p) {}
    DSAPub(DSA* pri); 
  };

  class DSA: public AbstractDSA {
  private:
    bint priKey;
  public:
    DSA(bool showoff = true); 
    DSA* setPri(bint pk) ;
    DSASign* sign(const bint msg); 
  };

  class DSASign {
  public:
    DSAPub* keyPub;
    bint msg, sig, r; 
    DSASign(DSAPub* keyPub, bint msg, bint sig, bint r):
      keyPub(keyPub),
      msg(msg),
      sig(sig),
      r(r) {}
    ~DSASign() {
      delete keyPub;
    }
    void showSign() {
      cout<<"\n Signature information: "<<endl;
      cout<<"Signer :"<<keyPub -> pubKey<<endl;
      cout<<"   Msg :"<<msg<<endl;
      cout<<"   Sig :"<<sig<<endl<<endl;
    }
    static void test();
  };

  // ------------------------------------

  void DSASign::test() {

    cout<<"\nSetting dsa signer 1:\n";
    DSA key1 = DSA();
    key1.setPri(bint("3894792874927"));

    cout<<"\nSetting dsa signer 2:\n";
    DSA key2 = DSA();
    key2.setPri(bint("8934798738927429"));

#ifdef DEBUG
    cout<<"key1.pub: "<<key1.pubKey<<endl;
    cout<<"key2.pub: "<<key2.pubKey<<endl;
#endif
    bint msg1 = bint("1234567890987654321");
    bint msg2 = bint("8888866666012345");
    DSASign *s1 = key1.sign(msg1);
    DSASign *s2 = key2.sign(msg1);
    DSASign *s3 = key2.sign(msg2);
    cout<<"\n key 1 signed msg1 while key 2 signed msg2 & msg3:\n\n";
    s1 -> showSign();
    s2 -> showSign();
    s3 -> showSign();
    cout<<"key 1 checks s1: "<<key1.check(s1) << endl;
    cout<<"key 1 checks s2: "<<key1.check(s2) << endl;
    cout<<"key 1 checks s3: "<<key1.check(s3) << endl;
    cout<<"key 2 checks s1: "<<key2.check(s1) << endl;
    cout<<"key 2 checks s2: "<<key2.check(s2) << endl;
    cout<<"key 2 checks s2: "<<key2.check(s3) << endl;
  }

  bool AbstractDSA::check(const DSASign* sObj) {
#ifdef PUBKEY_MATCH_BEFORE_SIGN_MATCH
    if (sObj -> keyPub-> pubKey!= pubKey) {
      cout<<"Blatant fabrication!\n";
      return false;
    }
    //#define DEBUG_checkPubkeysMatch
#ifdef DEBUG_checkPubkeysMatch
    else {
      cout<<"In sign: "<<sObj -> keyPub -> pubKey<<endl;
      cout<<"In key : "<<pubKey<<endl;
    }
#endif
#endif
    const bint &sig = sObj -> sig;
    const bint &r = sObj -> r;
    const bint &msg = sObj -> msg;
    const bint sig1 = modInverse(sig, q);
    //#define DEBUG_checkInverse
#ifdef DEBUG_checkInverse
    cout<<"check: muls inverse eqs 1: "<< (sig * sig1 % q) <<endl;
#endif
    const bint d1 = modpow(g, msg * sig1 % q, p);    // g^mv
    const bint d2 = modpow(pubKey, r * sig1 % q, p); // g^xrv
    const bint d = d1 * d2 % p;
#ifdef DEBUG_checkRsMeet
    cout<<"cheking d : "<<d<<endl;
    cout<<"  d mod q = "<<d % q<<endl;
    cout<<"while r = : "<<r<<endl;
#endif
    return (d1 * d2 % p) % q == r;
  }

  DSAPub::DSAPub(DSA *pri):
    DSAPub(pri -> pubKey, pri -> g, pri -> q, pri -> p)
  { }

  DSA::DSA(bool showoff) {
    q = goodPrime(bint("8232423489872349872394879089"));
    p = 2 * q ;
    bint k = 2;
    while (! isPrime(p + 1)) {
      p += 2 * q;
      k += 2;
    }
    p = p + 1;
    bint h = 6666;
    g = modpow(h, k, p);
    if (showoff) {
      cout<<"h = \t"<<h<<endl;
      cout<<"k = \t"<<k<<endl;
      cout<<"------- \n";
      cout<<"q = \t"<<q<<endl;
      cout<<"p = \t"<<p<<endl;
      cout<<"g = \t"<<g<<endl;
      cout<<"=======================\n\n";
    }
#ifdef DEBUG_relationsPQ
    {
      bint a = 191, b = 233;
      bint pub1 = modpow(g, a * b, p);
      bint pub2 = modpow(g, a * b % q, p);
      assert(pub1 == pub2);
      cout<<"\nPeriodic with p in  modpow of q"<<endl;
    }
#endif
  }

  DSA* DSA::setPri(bint pk) {
    priKey = pk % q;
    pubKey = modpow(g, priKey, p);
#ifdef DEBUG_privateSetReady
    cout<<"set: pri = "<<priKey<<" \t ";
    cout<<"pub = "<<pubKey<<endl;
#endif
    return this;
  }

  DSASign* DSA::sign(const bint msg) {
    bint k = modpow(msg, priKey, p) % q - 1;
    bint r = 0;
    while (r == 0) {
      k ++;
      r = modpow(g, k, p);
    } 
    r = r % q;

    const bint kInv = modInverse(k, q);
    const bint sign = ((msg + r * priKey % q) %q ) * kInv % q;
    DSAPub *pk = new DSAPub(this);
    DSASign *sObj = new DSASign (pk, msg, sign, r);
    //#define DEBUG_kInvIsRight
#ifdef DEBUG_kInvIsRight
    assert(modpow(sObj -> r, kInv, sObj -> keyPub -> p) == g);
    cout<<"r ^ kInv = g^(k*kInv) = g"<<endl;
#endif
#ifdef DEBUG_ifSignatureCanBeChecked
    const bint sInv = modInverse(sign, q);
    assert(modpow(g, (msg + priKey * r) * sInv % q, p) % q == r);
    assert(modpow(g, (msg ) * sInv % q, p) * modpow(g, priKey * r * sInv % q, p) % p % q == r);
    cout<<"So far so good...\n\n";
    // assert(modpow(g, (msg ) * sInv % q, p) * modpow(pubKey, r * sInv % q, p) % p % q == r);
    assert(modpow(g, priKey, p) == pubKey);
    /* Debug Conclusion: Public key should not modulo q ! */
#endif
    return sObj;
  }
}

namespace dh {
  class KeyPair{
  private:
    bint priKey;
    static bint epochPeriod ;
  public: 
    bint pubKey;
    bint base;
    bint modulus;
    KeyPair(bint pri, bint base, bint modu) {
      using namespace modcalcs;
      priKey = pri;
      this -> base = base;
      modulus = (modu);
      pubKey = modpow(base, pri, modulus);
    }
    KeyPair() {
      using namespace modcalcs;
      modulus = goodPrime(bint("3892492874398728934723248"));
      base = modulus / 3;
      for (; base < modulus - 8; base++) {
        if (modpow(base, (modulus - 1) / 2, modulus) != 1)
          break;
      }
      priKey = 1;
      pubKey = modpow(base, priKey, modulus);
    }
    KeyPair* copyFrom(const KeyPair& k1) {
      base = k1.base;
      modulus = k1.modulus;
      return this;
    }
    KeyPair* setPri(bint pri1) {
      using namespace modcalcs;
      this -> priKey = pri1;
      this -> pubKey = modpow(base, pri1, modulus);
      return this;
    }
    KeyPair* setPri(string pri1) {
      return setPri(bint(pri1));
    }
    bint exchange(KeyPair& bob) {
      using namespace modcalcs;
      return modpow(bob.pubKey, priKey, modulus);
    }
    bint crackLoop(const KeyPair* bob) {
      bint pri1 = 0, pub1 = 1;
      cout<<"Start cracking: "<<bob -> pubKey<<endl;
      for (; pub1 != bob -> pubKey; pri1++, pub1 = pub1 * base % modulus) {
        if (pri1 % epochPeriod == 1)
          cout<<pri1<< " failed to match ...    \r";
        if (pri1 >= modulus)
          return 0;
      }
      return pri1;
    }
    bint crackRecur(const KeyPair* bob) {
      if (pubKey == bob -> pubKey)
        return priKey;
      else if (priKey >= modulus)
        return priKey;
      else {
        if (priKey % epochPeriod == 1)
          cout<<priKey<< " still  doesn't work ...    \r";
        pubKey = pubKey * base % modulus;
        priKey++; 
        return crackRecur(bob);
      }
    }
    bint crack(const KeyPair* bob) {
      copyFrom(*bob);
      setPri(bint(0));
      cout<<"Cracking public key: "<<bob -> pubKey<<endl;
      return crackRecur(bob);
    }
    void showInfo(bool showpri = false) {
      cout<<"Modulus   :\t"<<modulus<<endl;
      cout<<"Base      :\t"<<base<<endl;
      cout<<"Pub key   :\t"<<pubKey<<endl;;
      if (showpri)
        cout<<"**Sec key!**  [\t"<<priKey<<" ]\n";
      cout<<endl;
    }

    static void test() {
      KeyPair a, b;
      a.setPri("28949287429834");
      b.copyFrom(a);
      b.setPri("512341234");
      cout<<"Alice's Key:\n";
      a.showInfo(true);
      cout<<"\nBob's Key:\n";
      b.showInfo(true);
      elapseTest(cout<<"Alice for Bob: "<<a.exchange(b)<<endl);
      elapseTest(cout<<"Bob for Alice: "<<b.exchange(a)<<endl);
      
      cout<<"\nAlice cracks Bob (Loops):\n\n";
      elapseTest(cout<<a.crackLoop(&b)<<"  success!   \n");

      cout<<"\nBob cracks Alice (Loops):\n\n";
      elapseTest(cout<<b.crackLoop(&a)<<"  success!   \n");
      // cout<<"\nAlice cracks Bob:\n\n";
      //elapseTest(cout<<a.crack(&b));
    }
  };
  bint KeyPair::epochPeriod = 1000000;
}


void modPowSpeed() {
  mpz_class m = mpz_class("23890423890273894273982891871"), prime;
  {
    using namespace modcalcs;
    cout<<"\n\n-- Recursion:\n";
    elapseTest(cout<<"result: "<<goodPrime(m)<<endl;);
    speeds::funcElapse([m] () {cout<<"lambda result: "<<goodPrime(m)<<endl;});
    // elapseTest(prime = goodPrime(m));
    // cout<<"result: "<<prime<<endl;
  }
  {
    using namespace modcalcs_iter;
    cout<<"\n\n-- Loops:\n";
    elapseTest(cout<<"result: "<<goodPrime(m)<<endl;);
    speeds::funcElapse([m] () {cout<<"lambda result: "<<goodPrime(m)<<endl;});
    //elapseTest(prime = goodPrime(m));
    //cout<<"result: "<<prime<<endl;
  }
}

int main() {
  using namespace toyDSA;
  DSASign::test();
  // modPowSpeed();
  // Bezouts::test();
  // dh::KeyPair::test();
  cout<<endl;;
  return 0;
}
