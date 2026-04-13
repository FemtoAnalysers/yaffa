// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME CreateDatasetDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "CreateDataset.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace Pythia8 {
   namespace ROOTDict {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *Pythia8_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("Pythia8", 0 /*version*/, "Pythia8/PythiaStdlib.h", 76,
                     ::ROOT::Internal::DefineBehavior((void*)nullptr,(void*)nullptr),
                     &Pythia8_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *Pythia8_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static TClass *Pythia8cLcLVec4_Dictionary();
   static void Pythia8cLcLVec4_TClassManip(TClass*);
   static void *new_Pythia8cLcLVec4(void *p = nullptr);
   static void *newArray_Pythia8cLcLVec4(Long_t size, void *p);
   static void delete_Pythia8cLcLVec4(void *p);
   static void deleteArray_Pythia8cLcLVec4(void *p);
   static void destruct_Pythia8cLcLVec4(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Vec4*)
   {
      ::Pythia8::Vec4 *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Vec4));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Vec4", "Pythia8/Basics.h", 32,
                  typeid(::Pythia8::Vec4), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLVec4_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Vec4) );
      instance.SetNew(&new_Pythia8cLcLVec4);
      instance.SetNewArray(&newArray_Pythia8cLcLVec4);
      instance.SetDelete(&delete_Pythia8cLcLVec4);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLVec4);
      instance.SetDestructor(&destruct_Pythia8cLcLVec4);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Vec4*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Vec4*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Vec4*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLVec4_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Vec4*>(nullptr))->GetClass();
      Pythia8cLcLVec4_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLVec4_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLRndmEngine_Dictionary();
   static void Pythia8cLcLRndmEngine_TClassManip(TClass*);
   static void *new_Pythia8cLcLRndmEngine(void *p = nullptr);
   static void *newArray_Pythia8cLcLRndmEngine(Long_t size, void *p);
   static void delete_Pythia8cLcLRndmEngine(void *p);
   static void deleteArray_Pythia8cLcLRndmEngine(void *p);
   static void destruct_Pythia8cLcLRndmEngine(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::RndmEngine*)
   {
      ::Pythia8::RndmEngine *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::RndmEngine));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::RndmEngine", "Pythia8/Basics.h", 352,
                  typeid(::Pythia8::RndmEngine), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLRndmEngine_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::RndmEngine) );
      instance.SetNew(&new_Pythia8cLcLRndmEngine);
      instance.SetNewArray(&newArray_Pythia8cLcLRndmEngine);
      instance.SetDelete(&delete_Pythia8cLcLRndmEngine);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLRndmEngine);
      instance.SetDestructor(&destruct_Pythia8cLcLRndmEngine);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::RndmEngine*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::RndmEngine*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::RndmEngine*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLRndmEngine_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::RndmEngine*>(nullptr))->GetClass();
      Pythia8cLcLRndmEngine_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLRndmEngine_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLRndm_Dictionary();
   static void Pythia8cLcLRndm_TClassManip(TClass*);
   static void *new_Pythia8cLcLRndm(void *p = nullptr);
   static void *newArray_Pythia8cLcLRndm(Long_t size, void *p);
   static void delete_Pythia8cLcLRndm(void *p);
   static void deleteArray_Pythia8cLcLRndm(void *p);
   static void destruct_Pythia8cLcLRndm(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Rndm*)
   {
      ::Pythia8::Rndm *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Rndm));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Rndm", "Pythia8/Basics.h", 385,
                  typeid(::Pythia8::Rndm), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLRndm_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Rndm) );
      instance.SetNew(&new_Pythia8cLcLRndm);
      instance.SetNewArray(&newArray_Pythia8cLcLRndm);
      instance.SetDelete(&delete_Pythia8cLcLRndm);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLRndm);
      instance.SetDestructor(&destruct_Pythia8cLcLRndm);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Rndm*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Rndm*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Rndm*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLRndm_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Rndm*>(nullptr))->GetClass();
      Pythia8cLcLRndm_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLRndm_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLInfo_Dictionary();
   static void Pythia8cLcLInfo_TClassManip(TClass*);
   static void *new_Pythia8cLcLInfo(void *p = nullptr);
   static void *newArray_Pythia8cLcLInfo(Long_t size, void *p);
   static void delete_Pythia8cLcLInfo(void *p);
   static void deleteArray_Pythia8cLcLInfo(void *p);
   static void destruct_Pythia8cLcLInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Info*)
   {
      ::Pythia8::Info *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Info));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Info", "Pythia8/Info.h", 45,
                  typeid(::Pythia8::Info), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLInfo_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Info) );
      instance.SetNew(&new_Pythia8cLcLInfo);
      instance.SetNewArray(&newArray_Pythia8cLcLInfo);
      instance.SetDelete(&delete_Pythia8cLcLInfo);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLInfo);
      instance.SetDestructor(&destruct_Pythia8cLcLInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Info*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Info*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Info*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Info*>(nullptr))->GetClass();
      Pythia8cLcLInfo_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLInfo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLSettings_Dictionary();
   static void Pythia8cLcLSettings_TClassManip(TClass*);
   static void *new_Pythia8cLcLSettings(void *p = nullptr);
   static void *newArray_Pythia8cLcLSettings(Long_t size, void *p);
   static void delete_Pythia8cLcLSettings(void *p);
   static void deleteArray_Pythia8cLcLSettings(void *p);
   static void destruct_Pythia8cLcLSettings(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Settings*)
   {
      ::Pythia8::Settings *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Settings));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Settings", "Pythia8/Settings.h", 195,
                  typeid(::Pythia8::Settings), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLSettings_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Settings) );
      instance.SetNew(&new_Pythia8cLcLSettings);
      instance.SetNewArray(&newArray_Pythia8cLcLSettings);
      instance.SetDelete(&delete_Pythia8cLcLSettings);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLSettings);
      instance.SetDestructor(&destruct_Pythia8cLcLSettings);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Settings*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Settings*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Settings*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLSettings_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Settings*>(nullptr))->GetClass();
      Pythia8cLcLSettings_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLSettings_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLAlphaStrong_Dictionary();
   static void Pythia8cLcLAlphaStrong_TClassManip(TClass*);
   static void *new_Pythia8cLcLAlphaStrong(void *p = nullptr);
   static void *newArray_Pythia8cLcLAlphaStrong(Long_t size, void *p);
   static void delete_Pythia8cLcLAlphaStrong(void *p);
   static void deleteArray_Pythia8cLcLAlphaStrong(void *p);
   static void destruct_Pythia8cLcLAlphaStrong(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::AlphaStrong*)
   {
      ::Pythia8::AlphaStrong *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::AlphaStrong));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::AlphaStrong", "Pythia8/StandardModel.h", 23,
                  typeid(::Pythia8::AlphaStrong), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLAlphaStrong_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::AlphaStrong) );
      instance.SetNew(&new_Pythia8cLcLAlphaStrong);
      instance.SetNewArray(&newArray_Pythia8cLcLAlphaStrong);
      instance.SetDelete(&delete_Pythia8cLcLAlphaStrong);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLAlphaStrong);
      instance.SetDestructor(&destruct_Pythia8cLcLAlphaStrong);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::AlphaStrong*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::AlphaStrong*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::AlphaStrong*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLAlphaStrong_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::AlphaStrong*>(nullptr))->GetClass();
      Pythia8cLcLAlphaStrong_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLAlphaStrong_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLAlphaEM_Dictionary();
   static void Pythia8cLcLAlphaEM_TClassManip(TClass*);
   static void *new_Pythia8cLcLAlphaEM(void *p = nullptr);
   static void *newArray_Pythia8cLcLAlphaEM(Long_t size, void *p);
   static void delete_Pythia8cLcLAlphaEM(void *p);
   static void deleteArray_Pythia8cLcLAlphaEM(void *p);
   static void destruct_Pythia8cLcLAlphaEM(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::AlphaEM*)
   {
      ::Pythia8::AlphaEM *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::AlphaEM));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::AlphaEM", "Pythia8/StandardModel.h", 106,
                  typeid(::Pythia8::AlphaEM), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLAlphaEM_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::AlphaEM) );
      instance.SetNew(&new_Pythia8cLcLAlphaEM);
      instance.SetNewArray(&newArray_Pythia8cLcLAlphaEM);
      instance.SetDelete(&delete_Pythia8cLcLAlphaEM);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLAlphaEM);
      instance.SetDestructor(&destruct_Pythia8cLcLAlphaEM);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::AlphaEM*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::AlphaEM*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::AlphaEM*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLAlphaEM_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::AlphaEM*>(nullptr))->GetClass();
      Pythia8cLcLAlphaEM_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLAlphaEM_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLCoupSM_Dictionary();
   static void Pythia8cLcLCoupSM_TClassManip(TClass*);
   static void *new_Pythia8cLcLCoupSM(void *p = nullptr);
   static void *newArray_Pythia8cLcLCoupSM(Long_t size, void *p);
   static void delete_Pythia8cLcLCoupSM(void *p);
   static void deleteArray_Pythia8cLcLCoupSM(void *p);
   static void destruct_Pythia8cLcLCoupSM(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::CoupSM*)
   {
      ::Pythia8::CoupSM *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::CoupSM));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::CoupSM", "Pythia8/StandardModel.h", 135,
                  typeid(::Pythia8::CoupSM), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLCoupSM_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::CoupSM) );
      instance.SetNew(&new_Pythia8cLcLCoupSM);
      instance.SetNewArray(&newArray_Pythia8cLcLCoupSM);
      instance.SetDelete(&delete_Pythia8cLcLCoupSM);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLCoupSM);
      instance.SetDestructor(&destruct_Pythia8cLcLCoupSM);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::CoupSM*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::CoupSM*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::CoupSM*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLCoupSM_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::CoupSM*>(nullptr))->GetClass();
      Pythia8cLcLCoupSM_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLCoupSM_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLDecayChannel_Dictionary();
   static void Pythia8cLcLDecayChannel_TClassManip(TClass*);
   static void *new_Pythia8cLcLDecayChannel(void *p = nullptr);
   static void *newArray_Pythia8cLcLDecayChannel(Long_t size, void *p);
   static void delete_Pythia8cLcLDecayChannel(void *p);
   static void deleteArray_Pythia8cLcLDecayChannel(void *p);
   static void destruct_Pythia8cLcLDecayChannel(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::DecayChannel*)
   {
      ::Pythia8::DecayChannel *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::DecayChannel));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::DecayChannel", "Pythia8/ParticleData.h", 35,
                  typeid(::Pythia8::DecayChannel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLDecayChannel_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::DecayChannel) );
      instance.SetNew(&new_Pythia8cLcLDecayChannel);
      instance.SetNewArray(&newArray_Pythia8cLcLDecayChannel);
      instance.SetDelete(&delete_Pythia8cLcLDecayChannel);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLDecayChannel);
      instance.SetDestructor(&destruct_Pythia8cLcLDecayChannel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::DecayChannel*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::DecayChannel*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::DecayChannel*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLDecayChannel_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::DecayChannel*>(nullptr))->GetClass();
      Pythia8cLcLDecayChannel_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLDecayChannel_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLParticleDataEntry_Dictionary();
   static void Pythia8cLcLParticleDataEntry_TClassManip(TClass*);
   static void *new_Pythia8cLcLParticleDataEntry(void *p = nullptr);
   static void *newArray_Pythia8cLcLParticleDataEntry(Long_t size, void *p);
   static void delete_Pythia8cLcLParticleDataEntry(void *p);
   static void deleteArray_Pythia8cLcLParticleDataEntry(void *p);
   static void destruct_Pythia8cLcLParticleDataEntry(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ParticleDataEntry*)
   {
      ::Pythia8::ParticleDataEntry *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ParticleDataEntry));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ParticleDataEntry", "Pythia8/ParticleData.h", 125,
                  typeid(::Pythia8::ParticleDataEntry), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLParticleDataEntry_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ParticleDataEntry) );
      instance.SetNew(&new_Pythia8cLcLParticleDataEntry);
      instance.SetNewArray(&newArray_Pythia8cLcLParticleDataEntry);
      instance.SetDelete(&delete_Pythia8cLcLParticleDataEntry);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLParticleDataEntry);
      instance.SetDestructor(&destruct_Pythia8cLcLParticleDataEntry);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ParticleDataEntry*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ParticleDataEntry*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ParticleDataEntry*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLParticleDataEntry_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ParticleDataEntry*>(nullptr))->GetClass();
      Pythia8cLcLParticleDataEntry_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLParticleDataEntry_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLParticleData_Dictionary();
   static void Pythia8cLcLParticleData_TClassManip(TClass*);
   static void *new_Pythia8cLcLParticleData(void *p = nullptr);
   static void *newArray_Pythia8cLcLParticleData(Long_t size, void *p);
   static void delete_Pythia8cLcLParticleData(void *p);
   static void deleteArray_Pythia8cLcLParticleData(void *p);
   static void destruct_Pythia8cLcLParticleData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ParticleData*)
   {
      ::Pythia8::ParticleData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ParticleData));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ParticleData", "Pythia8/ParticleData.h", 422,
                  typeid(::Pythia8::ParticleData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLParticleData_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ParticleData) );
      instance.SetNew(&new_Pythia8cLcLParticleData);
      instance.SetNewArray(&newArray_Pythia8cLcLParticleData);
      instance.SetDelete(&delete_Pythia8cLcLParticleData);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLParticleData);
      instance.SetDestructor(&destruct_Pythia8cLcLParticleData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ParticleData*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ParticleData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ParticleData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLParticleData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ParticleData*>(nullptr))->GetClass();
      Pythia8cLcLParticleData_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLParticleData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLParticle_Dictionary();
   static void Pythia8cLcLParticle_TClassManip(TClass*);
   static void *new_Pythia8cLcLParticle(void *p = nullptr);
   static void *newArray_Pythia8cLcLParticle(Long_t size, void *p);
   static void delete_Pythia8cLcLParticle(void *p);
   static void deleteArray_Pythia8cLcLParticle(void *p);
   static void destruct_Pythia8cLcLParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Particle*)
   {
      ::Pythia8::Particle *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Particle));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Particle", "Pythia8/Event.h", 32,
                  typeid(::Pythia8::Particle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLParticle_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Particle) );
      instance.SetNew(&new_Pythia8cLcLParticle);
      instance.SetNewArray(&newArray_Pythia8cLcLParticle);
      instance.SetDelete(&delete_Pythia8cLcLParticle);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLParticle);
      instance.SetDestructor(&destruct_Pythia8cLcLParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Particle*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Particle*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Particle*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLParticle_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Particle*>(nullptr))->GetClass();
      Pythia8cLcLParticle_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLParticle_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLJunction_Dictionary();
   static void Pythia8cLcLJunction_TClassManip(TClass*);
   static void *new_Pythia8cLcLJunction(void *p = nullptr);
   static void *newArray_Pythia8cLcLJunction(Long_t size, void *p);
   static void delete_Pythia8cLcLJunction(void *p);
   static void deleteArray_Pythia8cLcLJunction(void *p);
   static void destruct_Pythia8cLcLJunction(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Junction*)
   {
      ::Pythia8::Junction *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Junction));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Junction", "Pythia8/Event.h", 338,
                  typeid(::Pythia8::Junction), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLJunction_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Junction) );
      instance.SetNew(&new_Pythia8cLcLJunction);
      instance.SetNewArray(&newArray_Pythia8cLcLJunction);
      instance.SetDelete(&delete_Pythia8cLcLJunction);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLJunction);
      instance.SetDestructor(&destruct_Pythia8cLcLJunction);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Junction*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Junction*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Junction*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLJunction_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Junction*>(nullptr))->GetClass();
      Pythia8cLcLJunction_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLJunction_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLEvent_Dictionary();
   static void Pythia8cLcLEvent_TClassManip(TClass*);
   static void *new_Pythia8cLcLEvent(void *p = nullptr);
   static void *newArray_Pythia8cLcLEvent(Long_t size, void *p);
   static void delete_Pythia8cLcLEvent(void *p);
   static void deleteArray_Pythia8cLcLEvent(void *p);
   static void destruct_Pythia8cLcLEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Event*)
   {
      ::Pythia8::Event *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Event));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Event", "Pythia8/Event.h", 453,
                  typeid(::Pythia8::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLEvent_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Event) );
      instance.SetNew(&new_Pythia8cLcLEvent);
      instance.SetNewArray(&newArray_Pythia8cLcLEvent);
      instance.SetDelete(&delete_Pythia8cLcLEvent);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLEvent);
      instance.SetDestructor(&destruct_Pythia8cLcLEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Event*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Event*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Event*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Event*>(nullptr))->GetClass();
      Pythia8cLcLEvent_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceGeneric_Dictionary();
   static void Pythia8cLcLResonanceGeneric_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceGeneric(void *p);
   static void deleteArray_Pythia8cLcLResonanceGeneric(void *p);
   static void destruct_Pythia8cLcLResonanceGeneric(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceGeneric*)
   {
      ::Pythia8::ResonanceGeneric *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceGeneric));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceGeneric", "Pythia8/ResonanceWidths.h", 152,
                  typeid(::Pythia8::ResonanceGeneric), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceGeneric_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceGeneric) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceGeneric);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceGeneric);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceGeneric);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceGeneric*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceGeneric*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceGeneric*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceGeneric_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceGeneric*>(nullptr))->GetClass();
      Pythia8cLcLResonanceGeneric_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceGeneric_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceGmZ_Dictionary();
   static void Pythia8cLcLResonanceGmZ_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceGmZ(void *p);
   static void deleteArray_Pythia8cLcLResonanceGmZ(void *p);
   static void destruct_Pythia8cLcLResonanceGmZ(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceGmZ*)
   {
      ::Pythia8::ResonanceGmZ *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceGmZ));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceGmZ", "Pythia8/ResonanceWidths.h", 168,
                  typeid(::Pythia8::ResonanceGmZ), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceGmZ_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceGmZ) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceGmZ);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceGmZ);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceGmZ);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceGmZ*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceGmZ*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceGmZ*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceGmZ_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceGmZ*>(nullptr))->GetClass();
      Pythia8cLcLResonanceGmZ_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceGmZ_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceW_Dictionary();
   static void Pythia8cLcLResonanceW_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceW(void *p);
   static void deleteArray_Pythia8cLcLResonanceW(void *p);
   static void destruct_Pythia8cLcLResonanceW(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceW*)
   {
      ::Pythia8::ResonanceW *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceW));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceW", "Pythia8/ResonanceWidths.h", 197,
                  typeid(::Pythia8::ResonanceW), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceW_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceW) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceW);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceW);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceW);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceW*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceW*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceW*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceW_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceW*>(nullptr))->GetClass();
      Pythia8cLcLResonanceW_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceW_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceTop_Dictionary();
   static void Pythia8cLcLResonanceTop_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceTop(void *p);
   static void deleteArray_Pythia8cLcLResonanceTop(void *p);
   static void destruct_Pythia8cLcLResonanceTop(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceTop*)
   {
      ::Pythia8::ResonanceTop *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceTop));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceTop", "Pythia8/ResonanceWidths.h", 224,
                  typeid(::Pythia8::ResonanceTop), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceTop_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceTop) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceTop);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceTop);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceTop);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceTop*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceTop*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceTop*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceTop_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceTop*>(nullptr))->GetClass();
      Pythia8cLcLResonanceTop_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceTop_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceFour_Dictionary();
   static void Pythia8cLcLResonanceFour_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceFour(void *p);
   static void deleteArray_Pythia8cLcLResonanceFour(void *p);
   static void destruct_Pythia8cLcLResonanceFour(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceFour*)
   {
      ::Pythia8::ResonanceFour *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceFour));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceFour", "Pythia8/ResonanceWidths.h", 252,
                  typeid(::Pythia8::ResonanceFour), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceFour_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceFour) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceFour);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceFour);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceFour);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceFour*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceFour*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceFour*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceFour_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceFour*>(nullptr))->GetClass();
      Pythia8cLcLResonanceFour_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceFour_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceH_Dictionary();
   static void Pythia8cLcLResonanceH_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceH(void *p);
   static void deleteArray_Pythia8cLcLResonanceH(void *p);
   static void destruct_Pythia8cLcLResonanceH(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceH*)
   {
      ::Pythia8::ResonanceH *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceH));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceH", "Pythia8/ResonanceWidths.h", 280,
                  typeid(::Pythia8::ResonanceH), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceH_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceH) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceH);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceH);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceH);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceH*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceH*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceH*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceH_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceH*>(nullptr))->GetClass();
      Pythia8cLcLResonanceH_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceH_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceHchg_Dictionary();
   static void Pythia8cLcLResonanceHchg_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceHchg(void *p);
   static void deleteArray_Pythia8cLcLResonanceHchg(void *p);
   static void destruct_Pythia8cLcLResonanceHchg(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceHchg*)
   {
      ::Pythia8::ResonanceHchg *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceHchg));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceHchg", "Pythia8/ResonanceWidths.h", 333,
                  typeid(::Pythia8::ResonanceHchg), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceHchg_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceHchg) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceHchg);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceHchg);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceHchg);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceHchg*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceHchg*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceHchg*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceHchg_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceHchg*>(nullptr))->GetClass();
      Pythia8cLcLResonanceHchg_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceHchg_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceZprime_Dictionary();
   static void Pythia8cLcLResonanceZprime_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceZprime(void *p);
   static void deleteArray_Pythia8cLcLResonanceZprime(void *p);
   static void destruct_Pythia8cLcLResonanceZprime(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceZprime*)
   {
      ::Pythia8::ResonanceZprime *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceZprime));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceZprime", "Pythia8/ResonanceWidths.h", 362,
                  typeid(::Pythia8::ResonanceZprime), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceZprime_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceZprime) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceZprime);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceZprime);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceZprime);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceZprime*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceZprime*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceZprime*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceZprime_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceZprime*>(nullptr))->GetClass();
      Pythia8cLcLResonanceZprime_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceZprime_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceWprime_Dictionary();
   static void Pythia8cLcLResonanceWprime_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceWprime(void *p);
   static void deleteArray_Pythia8cLcLResonanceWprime(void *p);
   static void destruct_Pythia8cLcLResonanceWprime(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceWprime*)
   {
      ::Pythia8::ResonanceWprime *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceWprime));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceWprime", "Pythia8/ResonanceWidths.h", 395,
                  typeid(::Pythia8::ResonanceWprime), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceWprime_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceWprime) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceWprime);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceWprime);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceWprime);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceWprime*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceWprime*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceWprime*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceWprime_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceWprime*>(nullptr))->GetClass();
      Pythia8cLcLResonanceWprime_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceWprime_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceRhorizontal_Dictionary();
   static void Pythia8cLcLResonanceRhorizontal_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceRhorizontal(void *p);
   static void deleteArray_Pythia8cLcLResonanceRhorizontal(void *p);
   static void destruct_Pythia8cLcLResonanceRhorizontal(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceRhorizontal*)
   {
      ::Pythia8::ResonanceRhorizontal *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceRhorizontal));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceRhorizontal", "Pythia8/ResonanceWidths.h", 423,
                  typeid(::Pythia8::ResonanceRhorizontal), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceRhorizontal_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceRhorizontal) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceRhorizontal);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceRhorizontal);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceRhorizontal);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceRhorizontal*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceRhorizontal*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceRhorizontal*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceRhorizontal_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceRhorizontal*>(nullptr))->GetClass();
      Pythia8cLcLResonanceRhorizontal_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceRhorizontal_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceExcited_Dictionary();
   static void Pythia8cLcLResonanceExcited_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceExcited(void *p);
   static void deleteArray_Pythia8cLcLResonanceExcited(void *p);
   static void destruct_Pythia8cLcLResonanceExcited(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceExcited*)
   {
      ::Pythia8::ResonanceExcited *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceExcited));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceExcited", "Pythia8/ResonanceWidths.h", 450,
                  typeid(::Pythia8::ResonanceExcited), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceExcited_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceExcited) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceExcited);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceExcited);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceExcited);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceExcited*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceExcited*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceExcited*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceExcited_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceExcited*>(nullptr))->GetClass();
      Pythia8cLcLResonanceExcited_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceExcited_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceGraviton_Dictionary();
   static void Pythia8cLcLResonanceGraviton_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceGraviton(void *p);
   static void deleteArray_Pythia8cLcLResonanceGraviton(void *p);
   static void destruct_Pythia8cLcLResonanceGraviton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceGraviton*)
   {
      ::Pythia8::ResonanceGraviton *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceGraviton));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceGraviton", "Pythia8/ResonanceWidths.h", 478,
                  typeid(::Pythia8::ResonanceGraviton), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceGraviton_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceGraviton) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceGraviton);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceGraviton);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceGraviton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceGraviton*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceGraviton*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceGraviton*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceGraviton_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceGraviton*>(nullptr))->GetClass();
      Pythia8cLcLResonanceGraviton_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceGraviton_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceLeptoquark_Dictionary();
   static void Pythia8cLcLResonanceLeptoquark_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceLeptoquark(void *p);
   static void deleteArray_Pythia8cLcLResonanceLeptoquark(void *p);
   static void destruct_Pythia8cLcLResonanceLeptoquark(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceLeptoquark*)
   {
      ::Pythia8::ResonanceLeptoquark *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceLeptoquark));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceLeptoquark", "Pythia8/ResonanceWidths.h", 546,
                  typeid(::Pythia8::ResonanceLeptoquark), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceLeptoquark_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceLeptoquark) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceLeptoquark);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceLeptoquark);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceLeptoquark);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceLeptoquark*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceLeptoquark*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceLeptoquark*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceLeptoquark_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceLeptoquark*>(nullptr))->GetClass();
      Pythia8cLcLResonanceLeptoquark_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceLeptoquark_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceNuRight_Dictionary();
   static void Pythia8cLcLResonanceNuRight_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceNuRight(void *p);
   static void deleteArray_Pythia8cLcLResonanceNuRight(void *p);
   static void destruct_Pythia8cLcLResonanceNuRight(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceNuRight*)
   {
      ::Pythia8::ResonanceNuRight *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceNuRight));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceNuRight", "Pythia8/ResonanceWidths.h", 573,
                  typeid(::Pythia8::ResonanceNuRight), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceNuRight_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceNuRight) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceNuRight);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceNuRight);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceNuRight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceNuRight*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceNuRight*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceNuRight*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceNuRight_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceNuRight*>(nullptr))->GetClass();
      Pythia8cLcLResonanceNuRight_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceNuRight_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceZRight_Dictionary();
   static void Pythia8cLcLResonanceZRight_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceZRight(void *p);
   static void deleteArray_Pythia8cLcLResonanceZRight(void *p);
   static void destruct_Pythia8cLcLResonanceZRight(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceZRight*)
   {
      ::Pythia8::ResonanceZRight *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceZRight));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceZRight", "Pythia8/ResonanceWidths.h", 600,
                  typeid(::Pythia8::ResonanceZRight), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceZRight_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceZRight) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceZRight);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceZRight);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceZRight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceZRight*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceZRight*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceZRight*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceZRight_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceZRight*>(nullptr))->GetClass();
      Pythia8cLcLResonanceZRight_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceZRight_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceWRight_Dictionary();
   static void Pythia8cLcLResonanceWRight_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceWRight(void *p);
   static void deleteArray_Pythia8cLcLResonanceWRight(void *p);
   static void destruct_Pythia8cLcLResonanceWRight(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceWRight*)
   {
      ::Pythia8::ResonanceWRight *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceWRight));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceWRight", "Pythia8/ResonanceWidths.h", 627,
                  typeid(::Pythia8::ResonanceWRight), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceWRight_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceWRight) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceWRight);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceWRight);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceWRight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceWRight*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceWRight*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceWRight*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceWRight_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceWRight*>(nullptr))->GetClass();
      Pythia8cLcLResonanceWRight_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceWRight_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceHchgchgLeft_Dictionary();
   static void Pythia8cLcLResonanceHchgchgLeft_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceHchgchgLeft(void *p);
   static void deleteArray_Pythia8cLcLResonanceHchgchgLeft(void *p);
   static void destruct_Pythia8cLcLResonanceHchgchgLeft(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceHchgchgLeft*)
   {
      ::Pythia8::ResonanceHchgchgLeft *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceHchgchgLeft));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceHchgchgLeft", "Pythia8/ResonanceWidths.h", 654,
                  typeid(::Pythia8::ResonanceHchgchgLeft), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceHchgchgLeft_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceHchgchgLeft) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceHchgchgLeft);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceHchgchgLeft);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceHchgchgLeft);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceHchgchgLeft*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceHchgchgLeft*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceHchgchgLeft*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceHchgchgLeft_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceHchgchgLeft*>(nullptr))->GetClass();
      Pythia8cLcLResonanceHchgchgLeft_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceHchgchgLeft_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLResonanceHchgchgRight_Dictionary();
   static void Pythia8cLcLResonanceHchgchgRight_TClassManip(TClass*);
   static void delete_Pythia8cLcLResonanceHchgchgRight(void *p);
   static void deleteArray_Pythia8cLcLResonanceHchgchgRight(void *p);
   static void destruct_Pythia8cLcLResonanceHchgchgRight(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::ResonanceHchgchgRight*)
   {
      ::Pythia8::ResonanceHchgchgRight *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::ResonanceHchgchgRight));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::ResonanceHchgchgRight", "Pythia8/ResonanceWidths.h", 682,
                  typeid(::Pythia8::ResonanceHchgchgRight), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLResonanceHchgchgRight_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::ResonanceHchgchgRight) );
      instance.SetDelete(&delete_Pythia8cLcLResonanceHchgchgRight);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLResonanceHchgchgRight);
      instance.SetDestructor(&destruct_Pythia8cLcLResonanceHchgchgRight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::ResonanceHchgchgRight*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::ResonanceHchgchgRight*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceHchgchgRight*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLResonanceHchgchgRight_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::ResonanceHchgchgRight*>(nullptr))->GetClass();
      Pythia8cLcLResonanceHchgchgRight_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLResonanceHchgchgRight_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLInBeam_Dictionary();
   static void Pythia8cLcLInBeam_TClassManip(TClass*);
   static void *new_Pythia8cLcLInBeam(void *p = nullptr);
   static void *newArray_Pythia8cLcLInBeam(Long_t size, void *p);
   static void delete_Pythia8cLcLInBeam(void *p);
   static void deleteArray_Pythia8cLcLInBeam(void *p);
   static void destruct_Pythia8cLcLInBeam(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::InBeam*)
   {
      ::Pythia8::InBeam *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::InBeam));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::InBeam", "Pythia8/SigmaProcess.h", 51,
                  typeid(::Pythia8::InBeam), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLInBeam_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::InBeam) );
      instance.SetNew(&new_Pythia8cLcLInBeam);
      instance.SetNewArray(&newArray_Pythia8cLcLInBeam);
      instance.SetDelete(&delete_Pythia8cLcLInBeam);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLInBeam);
      instance.SetDestructor(&destruct_Pythia8cLcLInBeam);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::InBeam*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::InBeam*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::InBeam*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLInBeam_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::InBeam*>(nullptr))->GetClass();
      Pythia8cLcLInBeam_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLInBeam_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLInPair_Dictionary();
   static void Pythia8cLcLInPair_TClassManip(TClass*);
   static void *new_Pythia8cLcLInPair(void *p = nullptr);
   static void *newArray_Pythia8cLcLInPair(Long_t size, void *p);
   static void delete_Pythia8cLcLInPair(void *p);
   static void deleteArray_Pythia8cLcLInPair(void *p);
   static void destruct_Pythia8cLcLInPair(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::InPair*)
   {
      ::Pythia8::InPair *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::InPair));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::InPair", "Pythia8/SigmaProcess.h", 68,
                  typeid(::Pythia8::InPair), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLInPair_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::InPair) );
      instance.SetNew(&new_Pythia8cLcLInPair);
      instance.SetNewArray(&newArray_Pythia8cLcLInPair);
      instance.SetDelete(&delete_Pythia8cLcLInPair);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLInPair);
      instance.SetDestructor(&destruct_Pythia8cLcLInPair);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::InPair*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::InPair*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::InPair*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLInPair_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::InPair*>(nullptr))->GetClass();
      Pythia8cLcLInPair_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLInPair_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLSigmaProcess_Dictionary();
   static void Pythia8cLcLSigmaProcess_TClassManip(TClass*);
   static void delete_Pythia8cLcLSigmaProcess(void *p);
   static void deleteArray_Pythia8cLcLSigmaProcess(void *p);
   static void destruct_Pythia8cLcLSigmaProcess(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::SigmaProcess*)
   {
      ::Pythia8::SigmaProcess *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::SigmaProcess));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::SigmaProcess", "Pythia8/SigmaProcess.h", 86,
                  typeid(::Pythia8::SigmaProcess), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLSigmaProcess_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::SigmaProcess) );
      instance.SetDelete(&delete_Pythia8cLcLSigmaProcess);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLSigmaProcess);
      instance.SetDestructor(&destruct_Pythia8cLcLSigmaProcess);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::SigmaProcess*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::SigmaProcess*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::SigmaProcess*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLSigmaProcess_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::SigmaProcess*>(nullptr))->GetClass();
      Pythia8cLcLSigmaProcess_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLSigmaProcess_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *Pythia8cLcLPythia_Dictionary();
   static void Pythia8cLcLPythia_TClassManip(TClass*);
   static void *new_Pythia8cLcLPythia(void *p = nullptr);
   static void *newArray_Pythia8cLcLPythia(Long_t size, void *p);
   static void delete_Pythia8cLcLPythia(void *p);
   static void deleteArray_Pythia8cLcLPythia(void *p);
   static void destruct_Pythia8cLcLPythia(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Pythia8::Pythia*)
   {
      ::Pythia8::Pythia *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::Pythia8::Pythia));
      static ::ROOT::TGenericClassInfo 
         instance("Pythia8::Pythia", "Pythia8/Pythia.h", 71,
                  typeid(::Pythia8::Pythia), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &Pythia8cLcLPythia_Dictionary, isa_proxy, 4,
                  sizeof(::Pythia8::Pythia) );
      instance.SetNew(&new_Pythia8cLcLPythia);
      instance.SetNewArray(&newArray_Pythia8cLcLPythia);
      instance.SetDelete(&delete_Pythia8cLcLPythia);
      instance.SetDeleteArray(&deleteArray_Pythia8cLcLPythia);
      instance.SetDestructor(&destruct_Pythia8cLcLPythia);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Pythia8::Pythia*)
   {
      return GenerateInitInstanceLocal(static_cast<::Pythia8::Pythia*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Pythia8::Pythia*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *Pythia8cLcLPythia_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::Pythia8::Pythia*>(nullptr))->GetClass();
      Pythia8cLcLPythia_TClassManip(theClass);
   return theClass;
   }

   static void Pythia8cLcLPythia_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLVec4(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Vec4 : new ::Pythia8::Vec4;
   }
   static void *newArray_Pythia8cLcLVec4(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Vec4[nElements] : new ::Pythia8::Vec4[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLVec4(void *p) {
      delete (static_cast<::Pythia8::Vec4*>(p));
   }
   static void deleteArray_Pythia8cLcLVec4(void *p) {
      delete [] (static_cast<::Pythia8::Vec4*>(p));
   }
   static void destruct_Pythia8cLcLVec4(void *p) {
      typedef ::Pythia8::Vec4 current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Vec4

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLRndmEngine(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::RndmEngine : new ::Pythia8::RndmEngine;
   }
   static void *newArray_Pythia8cLcLRndmEngine(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::RndmEngine[nElements] : new ::Pythia8::RndmEngine[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLRndmEngine(void *p) {
      delete (static_cast<::Pythia8::RndmEngine*>(p));
   }
   static void deleteArray_Pythia8cLcLRndmEngine(void *p) {
      delete [] (static_cast<::Pythia8::RndmEngine*>(p));
   }
   static void destruct_Pythia8cLcLRndmEngine(void *p) {
      typedef ::Pythia8::RndmEngine current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::RndmEngine

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLRndm(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Rndm : new ::Pythia8::Rndm;
   }
   static void *newArray_Pythia8cLcLRndm(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Rndm[nElements] : new ::Pythia8::Rndm[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLRndm(void *p) {
      delete (static_cast<::Pythia8::Rndm*>(p));
   }
   static void deleteArray_Pythia8cLcLRndm(void *p) {
      delete [] (static_cast<::Pythia8::Rndm*>(p));
   }
   static void destruct_Pythia8cLcLRndm(void *p) {
      typedef ::Pythia8::Rndm current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Rndm

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLInfo(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Info : new ::Pythia8::Info;
   }
   static void *newArray_Pythia8cLcLInfo(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Info[nElements] : new ::Pythia8::Info[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLInfo(void *p) {
      delete (static_cast<::Pythia8::Info*>(p));
   }
   static void deleteArray_Pythia8cLcLInfo(void *p) {
      delete [] (static_cast<::Pythia8::Info*>(p));
   }
   static void destruct_Pythia8cLcLInfo(void *p) {
      typedef ::Pythia8::Info current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Info

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLSettings(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Settings : new ::Pythia8::Settings;
   }
   static void *newArray_Pythia8cLcLSettings(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Settings[nElements] : new ::Pythia8::Settings[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLSettings(void *p) {
      delete (static_cast<::Pythia8::Settings*>(p));
   }
   static void deleteArray_Pythia8cLcLSettings(void *p) {
      delete [] (static_cast<::Pythia8::Settings*>(p));
   }
   static void destruct_Pythia8cLcLSettings(void *p) {
      typedef ::Pythia8::Settings current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Settings

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLAlphaStrong(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::AlphaStrong : new ::Pythia8::AlphaStrong;
   }
   static void *newArray_Pythia8cLcLAlphaStrong(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::AlphaStrong[nElements] : new ::Pythia8::AlphaStrong[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLAlphaStrong(void *p) {
      delete (static_cast<::Pythia8::AlphaStrong*>(p));
   }
   static void deleteArray_Pythia8cLcLAlphaStrong(void *p) {
      delete [] (static_cast<::Pythia8::AlphaStrong*>(p));
   }
   static void destruct_Pythia8cLcLAlphaStrong(void *p) {
      typedef ::Pythia8::AlphaStrong current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::AlphaStrong

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLAlphaEM(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::AlphaEM : new ::Pythia8::AlphaEM;
   }
   static void *newArray_Pythia8cLcLAlphaEM(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::AlphaEM[nElements] : new ::Pythia8::AlphaEM[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLAlphaEM(void *p) {
      delete (static_cast<::Pythia8::AlphaEM*>(p));
   }
   static void deleteArray_Pythia8cLcLAlphaEM(void *p) {
      delete [] (static_cast<::Pythia8::AlphaEM*>(p));
   }
   static void destruct_Pythia8cLcLAlphaEM(void *p) {
      typedef ::Pythia8::AlphaEM current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::AlphaEM

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLCoupSM(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::CoupSM : new ::Pythia8::CoupSM;
   }
   static void *newArray_Pythia8cLcLCoupSM(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::CoupSM[nElements] : new ::Pythia8::CoupSM[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLCoupSM(void *p) {
      delete (static_cast<::Pythia8::CoupSM*>(p));
   }
   static void deleteArray_Pythia8cLcLCoupSM(void *p) {
      delete [] (static_cast<::Pythia8::CoupSM*>(p));
   }
   static void destruct_Pythia8cLcLCoupSM(void *p) {
      typedef ::Pythia8::CoupSM current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::CoupSM

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLDecayChannel(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::DecayChannel : new ::Pythia8::DecayChannel;
   }
   static void *newArray_Pythia8cLcLDecayChannel(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::DecayChannel[nElements] : new ::Pythia8::DecayChannel[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLDecayChannel(void *p) {
      delete (static_cast<::Pythia8::DecayChannel*>(p));
   }
   static void deleteArray_Pythia8cLcLDecayChannel(void *p) {
      delete [] (static_cast<::Pythia8::DecayChannel*>(p));
   }
   static void destruct_Pythia8cLcLDecayChannel(void *p) {
      typedef ::Pythia8::DecayChannel current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::DecayChannel

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLParticleDataEntry(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::ParticleDataEntry : new ::Pythia8::ParticleDataEntry;
   }
   static void *newArray_Pythia8cLcLParticleDataEntry(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::ParticleDataEntry[nElements] : new ::Pythia8::ParticleDataEntry[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLParticleDataEntry(void *p) {
      delete (static_cast<::Pythia8::ParticleDataEntry*>(p));
   }
   static void deleteArray_Pythia8cLcLParticleDataEntry(void *p) {
      delete [] (static_cast<::Pythia8::ParticleDataEntry*>(p));
   }
   static void destruct_Pythia8cLcLParticleDataEntry(void *p) {
      typedef ::Pythia8::ParticleDataEntry current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ParticleDataEntry

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLParticleData(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::ParticleData : new ::Pythia8::ParticleData;
   }
   static void *newArray_Pythia8cLcLParticleData(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::ParticleData[nElements] : new ::Pythia8::ParticleData[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLParticleData(void *p) {
      delete (static_cast<::Pythia8::ParticleData*>(p));
   }
   static void deleteArray_Pythia8cLcLParticleData(void *p) {
      delete [] (static_cast<::Pythia8::ParticleData*>(p));
   }
   static void destruct_Pythia8cLcLParticleData(void *p) {
      typedef ::Pythia8::ParticleData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ParticleData

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLParticle(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Particle : new ::Pythia8::Particle;
   }
   static void *newArray_Pythia8cLcLParticle(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Particle[nElements] : new ::Pythia8::Particle[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLParticle(void *p) {
      delete (static_cast<::Pythia8::Particle*>(p));
   }
   static void deleteArray_Pythia8cLcLParticle(void *p) {
      delete [] (static_cast<::Pythia8::Particle*>(p));
   }
   static void destruct_Pythia8cLcLParticle(void *p) {
      typedef ::Pythia8::Particle current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Particle

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLJunction(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Junction : new ::Pythia8::Junction;
   }
   static void *newArray_Pythia8cLcLJunction(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Junction[nElements] : new ::Pythia8::Junction[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLJunction(void *p) {
      delete (static_cast<::Pythia8::Junction*>(p));
   }
   static void deleteArray_Pythia8cLcLJunction(void *p) {
      delete [] (static_cast<::Pythia8::Junction*>(p));
   }
   static void destruct_Pythia8cLcLJunction(void *p) {
      typedef ::Pythia8::Junction current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Junction

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLEvent(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Event : new ::Pythia8::Event;
   }
   static void *newArray_Pythia8cLcLEvent(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Event[nElements] : new ::Pythia8::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLEvent(void *p) {
      delete (static_cast<::Pythia8::Event*>(p));
   }
   static void deleteArray_Pythia8cLcLEvent(void *p) {
      delete [] (static_cast<::Pythia8::Event*>(p));
   }
   static void destruct_Pythia8cLcLEvent(void *p) {
      typedef ::Pythia8::Event current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Event

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceGeneric(void *p) {
      delete (static_cast<::Pythia8::ResonanceGeneric*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceGeneric(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceGeneric*>(p));
   }
   static void destruct_Pythia8cLcLResonanceGeneric(void *p) {
      typedef ::Pythia8::ResonanceGeneric current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceGeneric

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceGmZ(void *p) {
      delete (static_cast<::Pythia8::ResonanceGmZ*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceGmZ(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceGmZ*>(p));
   }
   static void destruct_Pythia8cLcLResonanceGmZ(void *p) {
      typedef ::Pythia8::ResonanceGmZ current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceGmZ

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceW(void *p) {
      delete (static_cast<::Pythia8::ResonanceW*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceW(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceW*>(p));
   }
   static void destruct_Pythia8cLcLResonanceW(void *p) {
      typedef ::Pythia8::ResonanceW current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceW

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceTop(void *p) {
      delete (static_cast<::Pythia8::ResonanceTop*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceTop(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceTop*>(p));
   }
   static void destruct_Pythia8cLcLResonanceTop(void *p) {
      typedef ::Pythia8::ResonanceTop current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceTop

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceFour(void *p) {
      delete (static_cast<::Pythia8::ResonanceFour*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceFour(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceFour*>(p));
   }
   static void destruct_Pythia8cLcLResonanceFour(void *p) {
      typedef ::Pythia8::ResonanceFour current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceFour

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceH(void *p) {
      delete (static_cast<::Pythia8::ResonanceH*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceH(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceH*>(p));
   }
   static void destruct_Pythia8cLcLResonanceH(void *p) {
      typedef ::Pythia8::ResonanceH current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceH

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceHchg(void *p) {
      delete (static_cast<::Pythia8::ResonanceHchg*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceHchg(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceHchg*>(p));
   }
   static void destruct_Pythia8cLcLResonanceHchg(void *p) {
      typedef ::Pythia8::ResonanceHchg current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceHchg

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceZprime(void *p) {
      delete (static_cast<::Pythia8::ResonanceZprime*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceZprime(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceZprime*>(p));
   }
   static void destruct_Pythia8cLcLResonanceZprime(void *p) {
      typedef ::Pythia8::ResonanceZprime current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceZprime

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceWprime(void *p) {
      delete (static_cast<::Pythia8::ResonanceWprime*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceWprime(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceWprime*>(p));
   }
   static void destruct_Pythia8cLcLResonanceWprime(void *p) {
      typedef ::Pythia8::ResonanceWprime current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceWprime

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceRhorizontal(void *p) {
      delete (static_cast<::Pythia8::ResonanceRhorizontal*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceRhorizontal(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceRhorizontal*>(p));
   }
   static void destruct_Pythia8cLcLResonanceRhorizontal(void *p) {
      typedef ::Pythia8::ResonanceRhorizontal current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceRhorizontal

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceExcited(void *p) {
      delete (static_cast<::Pythia8::ResonanceExcited*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceExcited(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceExcited*>(p));
   }
   static void destruct_Pythia8cLcLResonanceExcited(void *p) {
      typedef ::Pythia8::ResonanceExcited current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceExcited

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceGraviton(void *p) {
      delete (static_cast<::Pythia8::ResonanceGraviton*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceGraviton(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceGraviton*>(p));
   }
   static void destruct_Pythia8cLcLResonanceGraviton(void *p) {
      typedef ::Pythia8::ResonanceGraviton current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceGraviton

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceLeptoquark(void *p) {
      delete (static_cast<::Pythia8::ResonanceLeptoquark*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceLeptoquark(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceLeptoquark*>(p));
   }
   static void destruct_Pythia8cLcLResonanceLeptoquark(void *p) {
      typedef ::Pythia8::ResonanceLeptoquark current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceLeptoquark

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceNuRight(void *p) {
      delete (static_cast<::Pythia8::ResonanceNuRight*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceNuRight(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceNuRight*>(p));
   }
   static void destruct_Pythia8cLcLResonanceNuRight(void *p) {
      typedef ::Pythia8::ResonanceNuRight current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceNuRight

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceZRight(void *p) {
      delete (static_cast<::Pythia8::ResonanceZRight*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceZRight(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceZRight*>(p));
   }
   static void destruct_Pythia8cLcLResonanceZRight(void *p) {
      typedef ::Pythia8::ResonanceZRight current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceZRight

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceWRight(void *p) {
      delete (static_cast<::Pythia8::ResonanceWRight*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceWRight(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceWRight*>(p));
   }
   static void destruct_Pythia8cLcLResonanceWRight(void *p) {
      typedef ::Pythia8::ResonanceWRight current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceWRight

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceHchgchgLeft(void *p) {
      delete (static_cast<::Pythia8::ResonanceHchgchgLeft*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceHchgchgLeft(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceHchgchgLeft*>(p));
   }
   static void destruct_Pythia8cLcLResonanceHchgchgLeft(void *p) {
      typedef ::Pythia8::ResonanceHchgchgLeft current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceHchgchgLeft

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLResonanceHchgchgRight(void *p) {
      delete (static_cast<::Pythia8::ResonanceHchgchgRight*>(p));
   }
   static void deleteArray_Pythia8cLcLResonanceHchgchgRight(void *p) {
      delete [] (static_cast<::Pythia8::ResonanceHchgchgRight*>(p));
   }
   static void destruct_Pythia8cLcLResonanceHchgchgRight(void *p) {
      typedef ::Pythia8::ResonanceHchgchgRight current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::ResonanceHchgchgRight

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLInBeam(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::InBeam : new ::Pythia8::InBeam;
   }
   static void *newArray_Pythia8cLcLInBeam(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::InBeam[nElements] : new ::Pythia8::InBeam[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLInBeam(void *p) {
      delete (static_cast<::Pythia8::InBeam*>(p));
   }
   static void deleteArray_Pythia8cLcLInBeam(void *p) {
      delete [] (static_cast<::Pythia8::InBeam*>(p));
   }
   static void destruct_Pythia8cLcLInBeam(void *p) {
      typedef ::Pythia8::InBeam current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::InBeam

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLInPair(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::InPair : new ::Pythia8::InPair;
   }
   static void *newArray_Pythia8cLcLInPair(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::InPair[nElements] : new ::Pythia8::InPair[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLInPair(void *p) {
      delete (static_cast<::Pythia8::InPair*>(p));
   }
   static void deleteArray_Pythia8cLcLInPair(void *p) {
      delete [] (static_cast<::Pythia8::InPair*>(p));
   }
   static void destruct_Pythia8cLcLInPair(void *p) {
      typedef ::Pythia8::InPair current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::InPair

namespace ROOT {
   // Wrapper around operator delete
   static void delete_Pythia8cLcLSigmaProcess(void *p) {
      delete (static_cast<::Pythia8::SigmaProcess*>(p));
   }
   static void deleteArray_Pythia8cLcLSigmaProcess(void *p) {
      delete [] (static_cast<::Pythia8::SigmaProcess*>(p));
   }
   static void destruct_Pythia8cLcLSigmaProcess(void *p) {
      typedef ::Pythia8::SigmaProcess current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::SigmaProcess

namespace ROOT {
   // Wrappers around operator new
   static void *new_Pythia8cLcLPythia(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Pythia : new ::Pythia8::Pythia;
   }
   static void *newArray_Pythia8cLcLPythia(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) ::Pythia8::Pythia[nElements] : new ::Pythia8::Pythia[nElements];
   }
   // Wrapper around operator delete
   static void delete_Pythia8cLcLPythia(void *p) {
      delete (static_cast<::Pythia8::Pythia*>(p));
   }
   static void deleteArray_Pythia8cLcLPythia(void *p) {
      delete [] (static_cast<::Pythia8::Pythia*>(p));
   }
   static void destruct_Pythia8cLcLPythia(void *p) {
      typedef ::Pythia8::Pythia current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Pythia8::Pythia

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = nullptr);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 425,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<string>","std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<string>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<string>*>(nullptr))->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete (static_cast<vector<string>*>(p));
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] (static_cast<vector<string>*>(p));
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR_Dictionary();
   static void vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(Long_t size, void *p);
   static void delete_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p);
   static void deleteArray_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p);
   static void destruct_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<shared_ptr<Pythia8::SigmaProcess> >*)
   {
      vector<shared_ptr<Pythia8::SigmaProcess> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<shared_ptr<Pythia8::SigmaProcess> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<shared_ptr<Pythia8::SigmaProcess> >", -2, "vector", 425,
                  typeid(vector<shared_ptr<Pythia8::SigmaProcess> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<shared_ptr<Pythia8::SigmaProcess> >) );
      instance.SetNew(&new_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR);
      instance.SetNewArray(&newArray_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR);
      instance.SetDelete(&delete_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR);
      instance.SetDestructor(&destruct_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<shared_ptr<Pythia8::SigmaProcess> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<shared_ptr<Pythia8::SigmaProcess> >","std::vector<std::shared_ptr<Pythia8::SigmaProcess>, std::allocator<std::shared_ptr<Pythia8::SigmaProcess> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<shared_ptr<Pythia8::SigmaProcess> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<shared_ptr<Pythia8::SigmaProcess> >*>(nullptr))->GetClass();
      vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<shared_ptr<Pythia8::SigmaProcess> > : new vector<shared_ptr<Pythia8::SigmaProcess> >;
   }
   static void *newArray_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<shared_ptr<Pythia8::SigmaProcess> >[nElements] : new vector<shared_ptr<Pythia8::SigmaProcess> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p) {
      delete (static_cast<vector<shared_ptr<Pythia8::SigmaProcess> >*>(p));
   }
   static void deleteArray_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p) {
      delete [] (static_cast<vector<shared_ptr<Pythia8::SigmaProcess> >*>(p));
   }
   static void destruct_vectorlEshared_ptrlEPythia8cLcLSigmaProcessgRsPgR(void *p) {
      typedef vector<shared_ptr<Pythia8::SigmaProcess> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<shared_ptr<Pythia8::SigmaProcess> >

namespace ROOT {
   static TClass *vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR_Dictionary();
   static void vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(Long_t size, void *p);
   static void delete_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p);
   static void deleteArray_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p);
   static void destruct_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<shared_ptr<Pythia8::ResonanceWidths> >*)
   {
      vector<shared_ptr<Pythia8::ResonanceWidths> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<shared_ptr<Pythia8::ResonanceWidths> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<shared_ptr<Pythia8::ResonanceWidths> >", -2, "vector", 425,
                  typeid(vector<shared_ptr<Pythia8::ResonanceWidths> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<shared_ptr<Pythia8::ResonanceWidths> >) );
      instance.SetNew(&new_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR);
      instance.SetNewArray(&newArray_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR);
      instance.SetDelete(&delete_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR);
      instance.SetDestructor(&destruct_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<shared_ptr<Pythia8::ResonanceWidths> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<shared_ptr<Pythia8::ResonanceWidths> >","std::vector<std::shared_ptr<Pythia8::ResonanceWidths>, std::allocator<std::shared_ptr<Pythia8::ResonanceWidths> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<shared_ptr<Pythia8::ResonanceWidths> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<shared_ptr<Pythia8::ResonanceWidths> >*>(nullptr))->GetClass();
      vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<shared_ptr<Pythia8::ResonanceWidths> > : new vector<shared_ptr<Pythia8::ResonanceWidths> >;
   }
   static void *newArray_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<shared_ptr<Pythia8::ResonanceWidths> >[nElements] : new vector<shared_ptr<Pythia8::ResonanceWidths> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p) {
      delete (static_cast<vector<shared_ptr<Pythia8::ResonanceWidths> >*>(p));
   }
   static void deleteArray_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p) {
      delete [] (static_cast<vector<shared_ptr<Pythia8::ResonanceWidths> >*>(p));
   }
   static void destruct_vectorlEshared_ptrlEPythia8cLcLResonanceWidthsgRsPgR(void *p) {
      typedef vector<shared_ptr<Pythia8::ResonanceWidths> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<shared_ptr<Pythia8::ResonanceWidths> >

namespace ROOT {
   static TClass *vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR_Dictionary();
   static void vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p = nullptr);
   static void *newArray_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(Long_t size, void *p);
   static void delete_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p);
   static void deleteArray_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p);
   static void destruct_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<shared_ptr<Pythia8::PhaseSpace> >*)
   {
      vector<shared_ptr<Pythia8::PhaseSpace> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<shared_ptr<Pythia8::PhaseSpace> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<shared_ptr<Pythia8::PhaseSpace> >", -2, "vector", 425,
                  typeid(vector<shared_ptr<Pythia8::PhaseSpace> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<shared_ptr<Pythia8::PhaseSpace> >) );
      instance.SetNew(&new_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR);
      instance.SetNewArray(&newArray_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR);
      instance.SetDelete(&delete_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR);
      instance.SetDestructor(&destruct_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<shared_ptr<Pythia8::PhaseSpace> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<shared_ptr<Pythia8::PhaseSpace> >","std::vector<std::shared_ptr<Pythia8::PhaseSpace>, std::allocator<std::shared_ptr<Pythia8::PhaseSpace> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<shared_ptr<Pythia8::PhaseSpace> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<shared_ptr<Pythia8::PhaseSpace> >*>(nullptr))->GetClass();
      vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<shared_ptr<Pythia8::PhaseSpace> > : new vector<shared_ptr<Pythia8::PhaseSpace> >;
   }
   static void *newArray_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<shared_ptr<Pythia8::PhaseSpace> >[nElements] : new vector<shared_ptr<Pythia8::PhaseSpace> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p) {
      delete (static_cast<vector<shared_ptr<Pythia8::PhaseSpace> >*>(p));
   }
   static void deleteArray_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p) {
      delete [] (static_cast<vector<shared_ptr<Pythia8::PhaseSpace> >*>(p));
   }
   static void destruct_vectorlEshared_ptrlEPythia8cLcLPhaseSpacegRsPgR(void *p) {
      typedef vector<shared_ptr<Pythia8::PhaseSpace> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<shared_ptr<Pythia8::PhaseSpace> >

namespace ROOT {
   static TClass *vectorlEpairlEintcOintgRsPgR_Dictionary();
   static void vectorlEpairlEintcOintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEintcOintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEpairlEintcOintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEintcOintgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEintcOintgRsPgR(void *p);
   static void destruct_vectorlEpairlEintcOintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<int,int> >*)
   {
      vector<pair<int,int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<int,int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<int,int> >", -2, "vector", 425,
                  typeid(vector<pair<int,int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEintcOintgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<int,int> >) );
      instance.SetNew(&new_vectorlEpairlEintcOintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEintcOintgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEintcOintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEintcOintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEintcOintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<int,int> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<pair<int,int> >","std::vector<std::pair<int, int>, std::allocator<std::pair<int, int> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<pair<int,int> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEintcOintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<pair<int,int> >*>(nullptr))->GetClass();
      vectorlEpairlEintcOintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEintcOintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEintcOintgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<int,int> > : new vector<pair<int,int> >;
   }
   static void *newArray_vectorlEpairlEintcOintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<pair<int,int> >[nElements] : new vector<pair<int,int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEintcOintgRsPgR(void *p) {
      delete (static_cast<vector<pair<int,int> >*>(p));
   }
   static void deleteArray_vectorlEpairlEintcOintgRsPgR(void *p) {
      delete [] (static_cast<vector<pair<int,int> >*>(p));
   }
   static void destruct_vectorlEpairlEintcOintgRsPgR(void *p) {
      typedef vector<pair<int,int> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<pair<int,int> >

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 425,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr))->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete (static_cast<vector<int>*>(p));
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] (static_cast<vector<int>*>(p));
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 425,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLVec4gR_Dictionary();
   static void vectorlEPythia8cLcLVec4gR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLVec4gR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLVec4gR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLVec4gR(void *p);
   static void deleteArray_vectorlEPythia8cLcLVec4gR(void *p);
   static void destruct_vectorlEPythia8cLcLVec4gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::Vec4>*)
   {
      vector<Pythia8::Vec4> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::Vec4>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::Vec4>", -2, "vector", 425,
                  typeid(vector<Pythia8::Vec4>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLVec4gR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::Vec4>) );
      instance.SetNew(&new_vectorlEPythia8cLcLVec4gR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLVec4gR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLVec4gR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLVec4gR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLVec4gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::Vec4> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::Vec4>","std::vector<Pythia8::Vec4, std::allocator<Pythia8::Vec4> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::Vec4>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLVec4gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::Vec4>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLVec4gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLVec4gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLVec4gR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::Vec4> : new vector<Pythia8::Vec4>;
   }
   static void *newArray_vectorlEPythia8cLcLVec4gR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::Vec4>[nElements] : new vector<Pythia8::Vec4>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLVec4gR(void *p) {
      delete (static_cast<vector<Pythia8::Vec4>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLVec4gR(void *p) {
      delete [] (static_cast<vector<Pythia8::Vec4>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLVec4gR(void *p) {
      typedef vector<Pythia8::Vec4> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::Vec4>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLStringBreaksgR_Dictionary();
   static void vectorlEPythia8cLcLStringBreaksgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLStringBreaksgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLStringBreaksgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLStringBreaksgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLStringBreaksgR(void *p);
   static void destruct_vectorlEPythia8cLcLStringBreaksgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::StringBreaks>*)
   {
      vector<Pythia8::StringBreaks> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::StringBreaks>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::StringBreaks>", -2, "vector", 425,
                  typeid(vector<Pythia8::StringBreaks>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLStringBreaksgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::StringBreaks>) );
      instance.SetNew(&new_vectorlEPythia8cLcLStringBreaksgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLStringBreaksgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLStringBreaksgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLStringBreaksgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLStringBreaksgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::StringBreaks> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::StringBreaks>","std::vector<Pythia8::StringBreaks, std::allocator<Pythia8::StringBreaks> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::StringBreaks>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLStringBreaksgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::StringBreaks>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLStringBreaksgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLStringBreaksgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLStringBreaksgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::StringBreaks> : new vector<Pythia8::StringBreaks>;
   }
   static void *newArray_vectorlEPythia8cLcLStringBreaksgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::StringBreaks>[nElements] : new vector<Pythia8::StringBreaks>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLStringBreaksgR(void *p) {
      delete (static_cast<vector<Pythia8::StringBreaks>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLStringBreaksgR(void *p) {
      delete [] (static_cast<vector<Pythia8::StringBreaks>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLStringBreaksgR(void *p) {
      typedef vector<Pythia8::StringBreaks> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::StringBreaks>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLPhysicsBasemUgR_Dictionary();
   static void vectorlEPythia8cLcLPhysicsBasemUgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLPhysicsBasemUgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLPhysicsBasemUgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLPhysicsBasemUgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLPhysicsBasemUgR(void *p);
   static void destruct_vectorlEPythia8cLcLPhysicsBasemUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::PhysicsBase*>*)
   {
      vector<Pythia8::PhysicsBase*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::PhysicsBase*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::PhysicsBase*>", -2, "vector", 425,
                  typeid(vector<Pythia8::PhysicsBase*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLPhysicsBasemUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::PhysicsBase*>) );
      instance.SetNew(&new_vectorlEPythia8cLcLPhysicsBasemUgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLPhysicsBasemUgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLPhysicsBasemUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLPhysicsBasemUgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLPhysicsBasemUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::PhysicsBase*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::PhysicsBase*>","std::vector<Pythia8::PhysicsBase*, std::allocator<Pythia8::PhysicsBase*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::PhysicsBase*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLPhysicsBasemUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::PhysicsBase*>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLPhysicsBasemUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLPhysicsBasemUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLPhysicsBasemUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::PhysicsBase*> : new vector<Pythia8::PhysicsBase*>;
   }
   static void *newArray_vectorlEPythia8cLcLPhysicsBasemUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::PhysicsBase*>[nElements] : new vector<Pythia8::PhysicsBase*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLPhysicsBasemUgR(void *p) {
      delete (static_cast<vector<Pythia8::PhysicsBase*>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLPhysicsBasemUgR(void *p) {
      delete [] (static_cast<vector<Pythia8::PhysicsBase*>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLPhysicsBasemUgR(void *p) {
      typedef vector<Pythia8::PhysicsBase*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::PhysicsBase*>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLParticlegR_Dictionary();
   static void vectorlEPythia8cLcLParticlegR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLParticlegR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLParticlegR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLParticlegR(void *p);
   static void deleteArray_vectorlEPythia8cLcLParticlegR(void *p);
   static void destruct_vectorlEPythia8cLcLParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::Particle>*)
   {
      vector<Pythia8::Particle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::Particle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::Particle>", -2, "vector", 425,
                  typeid(vector<Pythia8::Particle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLParticlegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::Particle>) );
      instance.SetNew(&new_vectorlEPythia8cLcLParticlegR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLParticlegR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLParticlegR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::Particle> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::Particle>","std::vector<Pythia8::Particle, std::allocator<Pythia8::Particle> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::Particle>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::Particle>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLParticlegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::Particle> : new vector<Pythia8::Particle>;
   }
   static void *newArray_vectorlEPythia8cLcLParticlegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::Particle>[nElements] : new vector<Pythia8::Particle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLParticlegR(void *p) {
      delete (static_cast<vector<Pythia8::Particle>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLParticlegR(void *p) {
      delete [] (static_cast<vector<Pythia8::Particle>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLParticlegR(void *p) {
      typedef vector<Pythia8::Particle> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::Particle>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLLHAgeneratorgR_Dictionary();
   static void vectorlEPythia8cLcLLHAgeneratorgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLLHAgeneratorgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLLHAgeneratorgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLLHAgeneratorgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLLHAgeneratorgR(void *p);
   static void destruct_vectorlEPythia8cLcLLHAgeneratorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::LHAgenerator>*)
   {
      vector<Pythia8::LHAgenerator> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::LHAgenerator>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::LHAgenerator>", -2, "vector", 425,
                  typeid(vector<Pythia8::LHAgenerator>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLLHAgeneratorgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::LHAgenerator>) );
      instance.SetNew(&new_vectorlEPythia8cLcLLHAgeneratorgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLLHAgeneratorgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLLHAgeneratorgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLLHAgeneratorgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLLHAgeneratorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::LHAgenerator> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::LHAgenerator>","std::vector<Pythia8::LHAgenerator, std::allocator<Pythia8::LHAgenerator> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::LHAgenerator>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLLHAgeneratorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::LHAgenerator>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLLHAgeneratorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLLHAgeneratorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLLHAgeneratorgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::LHAgenerator> : new vector<Pythia8::LHAgenerator>;
   }
   static void *newArray_vectorlEPythia8cLcLLHAgeneratorgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::LHAgenerator>[nElements] : new vector<Pythia8::LHAgenerator>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLLHAgeneratorgR(void *p) {
      delete (static_cast<vector<Pythia8::LHAgenerator>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLLHAgeneratorgR(void *p) {
      delete [] (static_cast<vector<Pythia8::LHAgenerator>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLLHAgeneratorgR(void *p) {
      typedef vector<Pythia8::LHAgenerator> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::LHAgenerator>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLJunctiongR_Dictionary();
   static void vectorlEPythia8cLcLJunctiongR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLJunctiongR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLJunctiongR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLJunctiongR(void *p);
   static void deleteArray_vectorlEPythia8cLcLJunctiongR(void *p);
   static void destruct_vectorlEPythia8cLcLJunctiongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::Junction>*)
   {
      vector<Pythia8::Junction> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::Junction>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::Junction>", -2, "vector", 425,
                  typeid(vector<Pythia8::Junction>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLJunctiongR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::Junction>) );
      instance.SetNew(&new_vectorlEPythia8cLcLJunctiongR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLJunctiongR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLJunctiongR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLJunctiongR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLJunctiongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::Junction> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::Junction>","std::vector<Pythia8::Junction, std::allocator<Pythia8::Junction> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::Junction>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLJunctiongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::Junction>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLJunctiongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLJunctiongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLJunctiongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::Junction> : new vector<Pythia8::Junction>;
   }
   static void *newArray_vectorlEPythia8cLcLJunctiongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::Junction>[nElements] : new vector<Pythia8::Junction>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLJunctiongR(void *p) {
      delete (static_cast<vector<Pythia8::Junction>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLJunctiongR(void *p) {
      delete [] (static_cast<vector<Pythia8::Junction>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLJunctiongR(void *p) {
      typedef vector<Pythia8::Junction> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::Junction>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLInPairgR_Dictionary();
   static void vectorlEPythia8cLcLInPairgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLInPairgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLInPairgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLInPairgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLInPairgR(void *p);
   static void destruct_vectorlEPythia8cLcLInPairgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::InPair>*)
   {
      vector<Pythia8::InPair> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::InPair>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::InPair>", -2, "vector", 425,
                  typeid(vector<Pythia8::InPair>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLInPairgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::InPair>) );
      instance.SetNew(&new_vectorlEPythia8cLcLInPairgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLInPairgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLInPairgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLInPairgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLInPairgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::InPair> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::InPair>","std::vector<Pythia8::InPair, std::allocator<Pythia8::InPair> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::InPair>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLInPairgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::InPair>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLInPairgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLInPairgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLInPairgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::InPair> : new vector<Pythia8::InPair>;
   }
   static void *newArray_vectorlEPythia8cLcLInPairgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::InPair>[nElements] : new vector<Pythia8::InPair>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLInPairgR(void *p) {
      delete (static_cast<vector<Pythia8::InPair>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLInPairgR(void *p) {
      delete [] (static_cast<vector<Pythia8::InPair>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLInPairgR(void *p) {
      typedef vector<Pythia8::InPair> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::InPair>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLInBeamgR_Dictionary();
   static void vectorlEPythia8cLcLInBeamgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLInBeamgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLInBeamgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLInBeamgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLInBeamgR(void *p);
   static void destruct_vectorlEPythia8cLcLInBeamgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::InBeam>*)
   {
      vector<Pythia8::InBeam> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::InBeam>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::InBeam>", -2, "vector", 425,
                  typeid(vector<Pythia8::InBeam>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLInBeamgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::InBeam>) );
      instance.SetNew(&new_vectorlEPythia8cLcLInBeamgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLInBeamgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLInBeamgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLInBeamgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLInBeamgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::InBeam> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::InBeam>","std::vector<Pythia8::InBeam, std::allocator<Pythia8::InBeam> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::InBeam>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLInBeamgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::InBeam>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLInBeamgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLInBeamgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLInBeamgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::InBeam> : new vector<Pythia8::InBeam>;
   }
   static void *newArray_vectorlEPythia8cLcLInBeamgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::InBeam>[nElements] : new vector<Pythia8::InBeam>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLInBeamgR(void *p) {
      delete (static_cast<vector<Pythia8::InBeam>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLInBeamgR(void *p) {
      delete [] (static_cast<vector<Pythia8::InBeam>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLInBeamgR(void *p) {
      typedef vector<Pythia8::InBeam> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::InBeam>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLHVcolsgR_Dictionary();
   static void vectorlEPythia8cLcLHVcolsgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLHVcolsgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLHVcolsgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLHVcolsgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLHVcolsgR(void *p);
   static void destruct_vectorlEPythia8cLcLHVcolsgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::HVcols>*)
   {
      vector<Pythia8::HVcols> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::HVcols>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::HVcols>", -2, "vector", 425,
                  typeid(vector<Pythia8::HVcols>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLHVcolsgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::HVcols>) );
      instance.SetNew(&new_vectorlEPythia8cLcLHVcolsgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLHVcolsgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLHVcolsgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLHVcolsgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLHVcolsgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::HVcols> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::HVcols>","std::vector<Pythia8::HVcols, std::allocator<Pythia8::HVcols> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::HVcols>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLHVcolsgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::HVcols>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLHVcolsgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLHVcolsgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLHVcolsgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::HVcols> : new vector<Pythia8::HVcols>;
   }
   static void *newArray_vectorlEPythia8cLcLHVcolsgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::HVcols>[nElements] : new vector<Pythia8::HVcols>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLHVcolsgR(void *p) {
      delete (static_cast<vector<Pythia8::HVcols>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLHVcolsgR(void *p) {
      delete [] (static_cast<vector<Pythia8::HVcols>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLHVcolsgR(void *p) {
      typedef vector<Pythia8::HVcols> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::HVcols>

namespace ROOT {
   static TClass *vectorlEPythia8cLcLDecayChannelgR_Dictionary();
   static void vectorlEPythia8cLcLDecayChannelgR_TClassManip(TClass*);
   static void *new_vectorlEPythia8cLcLDecayChannelgR(void *p = nullptr);
   static void *newArray_vectorlEPythia8cLcLDecayChannelgR(Long_t size, void *p);
   static void delete_vectorlEPythia8cLcLDecayChannelgR(void *p);
   static void deleteArray_vectorlEPythia8cLcLDecayChannelgR(void *p);
   static void destruct_vectorlEPythia8cLcLDecayChannelgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Pythia8::DecayChannel>*)
   {
      vector<Pythia8::DecayChannel> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Pythia8::DecayChannel>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Pythia8::DecayChannel>", -2, "vector", 425,
                  typeid(vector<Pythia8::DecayChannel>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEPythia8cLcLDecayChannelgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Pythia8::DecayChannel>) );
      instance.SetNew(&new_vectorlEPythia8cLcLDecayChannelgR);
      instance.SetNewArray(&newArray_vectorlEPythia8cLcLDecayChannelgR);
      instance.SetDelete(&delete_vectorlEPythia8cLcLDecayChannelgR);
      instance.SetDeleteArray(&deleteArray_vectorlEPythia8cLcLDecayChannelgR);
      instance.SetDestructor(&destruct_vectorlEPythia8cLcLDecayChannelgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Pythia8::DecayChannel> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<Pythia8::DecayChannel>","std::vector<Pythia8::DecayChannel, std::allocator<Pythia8::DecayChannel> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<Pythia8::DecayChannel>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEPythia8cLcLDecayChannelgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<Pythia8::DecayChannel>*>(nullptr))->GetClass();
      vectorlEPythia8cLcLDecayChannelgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEPythia8cLcLDecayChannelgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEPythia8cLcLDecayChannelgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::DecayChannel> : new vector<Pythia8::DecayChannel>;
   }
   static void *newArray_vectorlEPythia8cLcLDecayChannelgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<Pythia8::DecayChannel>[nElements] : new vector<Pythia8::DecayChannel>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEPythia8cLcLDecayChannelgR(void *p) {
      delete (static_cast<vector<Pythia8::DecayChannel>*>(p));
   }
   static void deleteArray_vectorlEPythia8cLcLDecayChannelgR(void *p) {
      delete [] (static_cast<vector<Pythia8::DecayChannel>*>(p));
   }
   static void destruct_vectorlEPythia8cLcLDecayChannelgR(void *p) {
      typedef vector<Pythia8::DecayChannel> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<Pythia8::DecayChannel>

namespace ROOT {
   static TClass *setlEstringgR_Dictionary();
   static void setlEstringgR_TClassManip(TClass*);
   static void *new_setlEstringgR(void *p = nullptr);
   static void *newArray_setlEstringgR(Long_t size, void *p);
   static void delete_setlEstringgR(void *p);
   static void deleteArray_setlEstringgR(void *p);
   static void destruct_setlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const set<string>*)
   {
      set<string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(set<string>));
      static ::ROOT::TGenericClassInfo 
         instance("set<string>", -2, "set", 96,
                  typeid(set<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &setlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(set<string>) );
      instance.SetNew(&new_setlEstringgR);
      instance.SetNewArray(&newArray_setlEstringgR);
      instance.SetDelete(&delete_setlEstringgR);
      instance.SetDeleteArray(&deleteArray_setlEstringgR);
      instance.SetDestructor(&destruct_setlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Insert< set<string> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("set<string>","std::set<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const set<string>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *setlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const set<string>*>(nullptr))->GetClass();
      setlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void setlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_setlEstringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) set<string> : new set<string>;
   }
   static void *newArray_setlEstringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) set<string>[nElements] : new set<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_setlEstringgR(void *p) {
      delete (static_cast<set<string>*>(p));
   }
   static void deleteArray_setlEstringgR(void *p) {
      delete [] (static_cast<set<string>*>(p));
   }
   static void destruct_setlEstringgR(void *p) {
      typedef set<string> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class set<string>

namespace ROOT {
   static TClass *maplEstringcOstringgR_Dictionary();
   static void maplEstringcOstringgR_TClassManip(TClass*);
   static void *new_maplEstringcOstringgR(void *p = nullptr);
   static void *newArray_maplEstringcOstringgR(Long_t size, void *p);
   static void delete_maplEstringcOstringgR(void *p);
   static void deleteArray_maplEstringcOstringgR(void *p);
   static void destruct_maplEstringcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,string>*)
   {
      map<string,string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,string>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,string>", -2, "map", 102,
                  typeid(map<string,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOstringgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,string>) );
      instance.SetNew(&new_maplEstringcOstringgR);
      instance.SetNewArray(&newArray_maplEstringcOstringgR);
      instance.SetDelete(&delete_maplEstringcOstringgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOstringgR);
      instance.SetDestructor(&destruct_maplEstringcOstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,string> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,string>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,string>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,string>*>(nullptr))->GetClass();
      maplEstringcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOstringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,string> : new map<string,string>;
   }
   static void *newArray_maplEstringcOstringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,string>[nElements] : new map<string,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOstringgR(void *p) {
      delete (static_cast<map<string,string>*>(p));
   }
   static void deleteArray_maplEstringcOstringgR(void *p) {
      delete [] (static_cast<map<string,string>*>(p));
   }
   static void destruct_maplEstringcOstringgR(void *p) {
      typedef map<string,string> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,string>

namespace ROOT {
   static TClass *maplEstringcOdoublegR_Dictionary();
   static void maplEstringcOdoublegR_TClassManip(TClass*);
   static void *new_maplEstringcOdoublegR(void *p = nullptr);
   static void *newArray_maplEstringcOdoublegR(Long_t size, void *p);
   static void delete_maplEstringcOdoublegR(void *p);
   static void deleteArray_maplEstringcOdoublegR(void *p);
   static void destruct_maplEstringcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,double>*)
   {
      map<string,double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,double>", -2, "map", 102,
                  typeid(map<string,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,double>) );
      instance.SetNew(&new_maplEstringcOdoublegR);
      instance.SetNewArray(&newArray_maplEstringcOdoublegR);
      instance.SetDelete(&delete_maplEstringcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOdoublegR);
      instance.SetDestructor(&destruct_maplEstringcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,double>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, double, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, double> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,double>*>(nullptr))->GetClass();
      maplEstringcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,double> : new map<string,double>;
   }
   static void *newArray_maplEstringcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,double>[nElements] : new map<string,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOdoublegR(void *p) {
      delete (static_cast<map<string,double>*>(p));
   }
   static void deleteArray_maplEstringcOdoublegR(void *p) {
      delete [] (static_cast<map<string,double>*>(p));
   }
   static void destruct_maplEstringcOdoublegR(void *p) {
      typedef map<string,double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,double>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLWordgR_Dictionary();
   static void maplEstringcOPythia8cLcLWordgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLWordgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLWordgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLWordgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLWordgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLWordgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::Word>*)
   {
      map<string,Pythia8::Word> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::Word>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::Word>", -2, "map", 102,
                  typeid(map<string,Pythia8::Word>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLWordgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::Word>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLWordgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLWordgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLWordgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLWordgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLWordgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::Word> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::Word>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::Word, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::Word> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Word>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLWordgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Word>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLWordgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLWordgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLWordgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Word> : new map<string,Pythia8::Word>;
   }
   static void *newArray_maplEstringcOPythia8cLcLWordgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Word>[nElements] : new map<string,Pythia8::Word>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLWordgR(void *p) {
      delete (static_cast<map<string,Pythia8::Word>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLWordgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::Word>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLWordgR(void *p) {
      typedef map<string,Pythia8::Word> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::Word>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLWVecgR_Dictionary();
   static void maplEstringcOPythia8cLcLWVecgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLWVecgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLWVecgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLWVecgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLWVecgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLWVecgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::WVec>*)
   {
      map<string,Pythia8::WVec> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::WVec>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::WVec>", -2, "map", 102,
                  typeid(map<string,Pythia8::WVec>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLWVecgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::WVec>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLWVecgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLWVecgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLWVecgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLWVecgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLWVecgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::WVec> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::WVec>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::WVec, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::WVec> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::WVec>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLWVecgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::WVec>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLWVecgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLWVecgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLWVecgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::WVec> : new map<string,Pythia8::WVec>;
   }
   static void *newArray_maplEstringcOPythia8cLcLWVecgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::WVec>[nElements] : new map<string,Pythia8::WVec>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLWVecgR(void *p) {
      delete (static_cast<map<string,Pythia8::WVec>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLWVecgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::WVec>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLWVecgR(void *p) {
      typedef map<string,Pythia8::WVec> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::WVec>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLParmgR_Dictionary();
   static void maplEstringcOPythia8cLcLParmgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLParmgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLParmgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLParmgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLParmgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLParmgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::Parm>*)
   {
      map<string,Pythia8::Parm> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::Parm>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::Parm>", -2, "map", 102,
                  typeid(map<string,Pythia8::Parm>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLParmgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::Parm>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLParmgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLParmgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLParmgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLParmgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLParmgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::Parm> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::Parm>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::Parm, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::Parm> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Parm>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLParmgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Parm>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLParmgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLParmgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLParmgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Parm> : new map<string,Pythia8::Parm>;
   }
   static void *newArray_maplEstringcOPythia8cLcLParmgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Parm>[nElements] : new map<string,Pythia8::Parm>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLParmgR(void *p) {
      delete (static_cast<map<string,Pythia8::Parm>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLParmgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::Parm>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLParmgR(void *p) {
      typedef map<string,Pythia8::Parm> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::Parm>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLPVecgR_Dictionary();
   static void maplEstringcOPythia8cLcLPVecgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLPVecgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLPVecgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLPVecgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLPVecgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLPVecgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::PVec>*)
   {
      map<string,Pythia8::PVec> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::PVec>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::PVec>", -2, "map", 102,
                  typeid(map<string,Pythia8::PVec>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLPVecgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::PVec>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLPVecgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLPVecgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLPVecgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLPVecgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLPVecgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::PVec> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::PVec>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::PVec, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::PVec> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::PVec>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLPVecgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::PVec>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLPVecgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLPVecgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLPVecgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::PVec> : new map<string,Pythia8::PVec>;
   }
   static void *newArray_maplEstringcOPythia8cLcLPVecgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::PVec>[nElements] : new map<string,Pythia8::PVec>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLPVecgR(void *p) {
      delete (static_cast<map<string,Pythia8::PVec>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLPVecgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::PVec>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLPVecgR(void *p) {
      typedef map<string,Pythia8::PVec> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::PVec>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLModegR_Dictionary();
   static void maplEstringcOPythia8cLcLModegR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLModegR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLModegR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLModegR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLModegR(void *p);
   static void destruct_maplEstringcOPythia8cLcLModegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::Mode>*)
   {
      map<string,Pythia8::Mode> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::Mode>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::Mode>", -2, "map", 102,
                  typeid(map<string,Pythia8::Mode>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLModegR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::Mode>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLModegR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLModegR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLModegR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLModegR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLModegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::Mode> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::Mode>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::Mode, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::Mode> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Mode>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLModegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Mode>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLModegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLModegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLModegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Mode> : new map<string,Pythia8::Mode>;
   }
   static void *newArray_maplEstringcOPythia8cLcLModegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Mode>[nElements] : new map<string,Pythia8::Mode>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLModegR(void *p) {
      delete (static_cast<map<string,Pythia8::Mode>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLModegR(void *p) {
      delete [] (static_cast<map<string,Pythia8::Mode>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLModegR(void *p) {
      typedef map<string,Pythia8::Mode> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::Mode>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLMVecgR_Dictionary();
   static void maplEstringcOPythia8cLcLMVecgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLMVecgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLMVecgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLMVecgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLMVecgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLMVecgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::MVec>*)
   {
      map<string,Pythia8::MVec> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::MVec>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::MVec>", -2, "map", 102,
                  typeid(map<string,Pythia8::MVec>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLMVecgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::MVec>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLMVecgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLMVecgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLMVecgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLMVecgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLMVecgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::MVec> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::MVec>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::MVec, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::MVec> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::MVec>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLMVecgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::MVec>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLMVecgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLMVecgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLMVecgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::MVec> : new map<string,Pythia8::MVec>;
   }
   static void *newArray_maplEstringcOPythia8cLcLMVecgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::MVec>[nElements] : new map<string,Pythia8::MVec>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLMVecgR(void *p) {
      delete (static_cast<map<string,Pythia8::MVec>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLMVecgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::MVec>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLMVecgR(void *p) {
      typedef map<string,Pythia8::MVec> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::MVec>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLLHAweightgroupgR_Dictionary();
   static void maplEstringcOPythia8cLcLLHAweightgroupgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLLHAweightgroupgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::LHAweightgroup>*)
   {
      map<string,Pythia8::LHAweightgroup> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::LHAweightgroup>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::LHAweightgroup>", -2, "map", 102,
                  typeid(map<string,Pythia8::LHAweightgroup>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLLHAweightgroupgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::LHAweightgroup>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLLHAweightgroupgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLLHAweightgroupgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLLHAweightgroupgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLLHAweightgroupgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLLHAweightgroupgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::LHAweightgroup> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::LHAweightgroup>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::LHAweightgroup, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::LHAweightgroup> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::LHAweightgroup>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLLHAweightgroupgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::LHAweightgroup>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLLHAweightgroupgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLLHAweightgroupgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::LHAweightgroup> : new map<string,Pythia8::LHAweightgroup>;
   }
   static void *newArray_maplEstringcOPythia8cLcLLHAweightgroupgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::LHAweightgroup>[nElements] : new map<string,Pythia8::LHAweightgroup>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p) {
      delete (static_cast<map<string,Pythia8::LHAweightgroup>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::LHAweightgroup>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLLHAweightgroupgR(void *p) {
      typedef map<string,Pythia8::LHAweightgroup> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::LHAweightgroup>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLLHAweightgR_Dictionary();
   static void maplEstringcOPythia8cLcLLHAweightgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLLHAweightgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLLHAweightgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLLHAweightgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLLHAweightgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLLHAweightgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::LHAweight>*)
   {
      map<string,Pythia8::LHAweight> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::LHAweight>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::LHAweight>", -2, "map", 102,
                  typeid(map<string,Pythia8::LHAweight>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLLHAweightgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::LHAweight>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLLHAweightgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLLHAweightgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLLHAweightgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLLHAweightgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLLHAweightgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::LHAweight> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::LHAweight>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::LHAweight, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::LHAweight> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::LHAweight>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLLHAweightgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::LHAweight>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLLHAweightgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLLHAweightgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLLHAweightgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::LHAweight> : new map<string,Pythia8::LHAweight>;
   }
   static void *newArray_maplEstringcOPythia8cLcLLHAweightgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::LHAweight>[nElements] : new map<string,Pythia8::LHAweight>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLLHAweightgR(void *p) {
      delete (static_cast<map<string,Pythia8::LHAweight>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLLHAweightgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::LHAweight>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLLHAweightgR(void *p) {
      typedef map<string,Pythia8::LHAweight> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::LHAweight>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLFlaggR_Dictionary();
   static void maplEstringcOPythia8cLcLFlaggR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLFlaggR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLFlaggR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLFlaggR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLFlaggR(void *p);
   static void destruct_maplEstringcOPythia8cLcLFlaggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::Flag>*)
   {
      map<string,Pythia8::Flag> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::Flag>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::Flag>", -2, "map", 102,
                  typeid(map<string,Pythia8::Flag>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLFlaggR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::Flag>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLFlaggR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLFlaggR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLFlaggR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLFlaggR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLFlaggR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::Flag> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::Flag>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::Flag, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::Flag> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Flag>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLFlaggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::Flag>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLFlaggR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLFlaggR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLFlaggR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Flag> : new map<string,Pythia8::Flag>;
   }
   static void *newArray_maplEstringcOPythia8cLcLFlaggR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::Flag>[nElements] : new map<string,Pythia8::Flag>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLFlaggR(void *p) {
      delete (static_cast<map<string,Pythia8::Flag>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLFlaggR(void *p) {
      delete [] (static_cast<map<string,Pythia8::Flag>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLFlaggR(void *p) {
      typedef map<string,Pythia8::Flag> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::Flag>

namespace ROOT {
   static TClass *maplEstringcOPythia8cLcLFVecgR_Dictionary();
   static void maplEstringcOPythia8cLcLFVecgR_TClassManip(TClass*);
   static void *new_maplEstringcOPythia8cLcLFVecgR(void *p = nullptr);
   static void *newArray_maplEstringcOPythia8cLcLFVecgR(Long_t size, void *p);
   static void delete_maplEstringcOPythia8cLcLFVecgR(void *p);
   static void deleteArray_maplEstringcOPythia8cLcLFVecgR(void *p);
   static void destruct_maplEstringcOPythia8cLcLFVecgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,Pythia8::FVec>*)
   {
      map<string,Pythia8::FVec> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,Pythia8::FVec>));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,Pythia8::FVec>", -2, "map", 102,
                  typeid(map<string,Pythia8::FVec>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOPythia8cLcLFVecgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,Pythia8::FVec>) );
      instance.SetNew(&new_maplEstringcOPythia8cLcLFVecgR);
      instance.SetNewArray(&newArray_maplEstringcOPythia8cLcLFVecgR);
      instance.SetDelete(&delete_maplEstringcOPythia8cLcLFVecgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOPythia8cLcLFVecgR);
      instance.SetDestructor(&destruct_maplEstringcOPythia8cLcLFVecgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,Pythia8::FVec> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<string,Pythia8::FVec>","std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Pythia8::FVec, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, Pythia8::FVec> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::FVec>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOPythia8cLcLFVecgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<string,Pythia8::FVec>*>(nullptr))->GetClass();
      maplEstringcOPythia8cLcLFVecgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOPythia8cLcLFVecgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOPythia8cLcLFVecgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::FVec> : new map<string,Pythia8::FVec>;
   }
   static void *newArray_maplEstringcOPythia8cLcLFVecgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<string,Pythia8::FVec>[nElements] : new map<string,Pythia8::FVec>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOPythia8cLcLFVecgR(void *p) {
      delete (static_cast<map<string,Pythia8::FVec>*>(p));
   }
   static void deleteArray_maplEstringcOPythia8cLcLFVecgR(void *p) {
      delete [] (static_cast<map<string,Pythia8::FVec>*>(p));
   }
   static void destruct_maplEstringcOPythia8cLcLFVecgR(void *p) {
      typedef map<string,Pythia8::FVec> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<string,Pythia8::FVec>

namespace ROOT {
   static TClass *maplEintcOvectorlEstringgRsPgR_Dictionary();
   static void maplEintcOvectorlEstringgRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOvectorlEstringgRsPgR(void *p = nullptr);
   static void *newArray_maplEintcOvectorlEstringgRsPgR(Long_t size, void *p);
   static void delete_maplEintcOvectorlEstringgRsPgR(void *p);
   static void deleteArray_maplEintcOvectorlEstringgRsPgR(void *p);
   static void destruct_maplEintcOvectorlEstringgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,vector<string> >*)
   {
      map<int,vector<string> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,vector<string> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,vector<string> >", -2, "map", 102,
                  typeid(map<int,vector<string> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOvectorlEstringgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,vector<string> >) );
      instance.SetNew(&new_maplEintcOvectorlEstringgRsPgR);
      instance.SetNewArray(&newArray_maplEintcOvectorlEstringgRsPgR);
      instance.SetDelete(&delete_maplEintcOvectorlEstringgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOvectorlEstringgRsPgR);
      instance.SetDestructor(&destruct_maplEintcOvectorlEstringgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,vector<string> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<int,vector<string> >","std::map<int, std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<int>, std::allocator<std::pair<int const, std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<int,vector<string> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOvectorlEstringgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<int,vector<string> >*>(nullptr))->GetClass();
      maplEintcOvectorlEstringgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOvectorlEstringgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOvectorlEstringgRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,vector<string> > : new map<int,vector<string> >;
   }
   static void *newArray_maplEintcOvectorlEstringgRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,vector<string> >[nElements] : new map<int,vector<string> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOvectorlEstringgRsPgR(void *p) {
      delete (static_cast<map<int,vector<string> >*>(p));
   }
   static void deleteArray_maplEintcOvectorlEstringgRsPgR(void *p) {
      delete [] (static_cast<map<int,vector<string> >*>(p));
   }
   static void destruct_maplEintcOvectorlEstringgRsPgR(void *p) {
      typedef map<int,vector<string> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<int,vector<string> >

namespace ROOT {
   static TClass *maplEintcOstringgR_Dictionary();
   static void maplEintcOstringgR_TClassManip(TClass*);
   static void *new_maplEintcOstringgR(void *p = nullptr);
   static void *newArray_maplEintcOstringgR(Long_t size, void *p);
   static void delete_maplEintcOstringgR(void *p);
   static void deleteArray_maplEintcOstringgR(void *p);
   static void destruct_maplEintcOstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,string>*)
   {
      map<int,string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,string>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,string>", -2, "map", 102,
                  typeid(map<int,string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOstringgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,string>) );
      instance.SetNew(&new_maplEintcOstringgR);
      instance.SetNewArray(&newArray_maplEintcOstringgR);
      instance.SetDelete(&delete_maplEintcOstringgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOstringgR);
      instance.SetDestructor(&destruct_maplEintcOstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,string> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<int,string>","std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<int,string>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<int,string>*>(nullptr))->GetClass();
      maplEintcOstringgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOstringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,string> : new map<int,string>;
   }
   static void *newArray_maplEintcOstringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,string>[nElements] : new map<int,string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOstringgR(void *p) {
      delete (static_cast<map<int,string>*>(p));
   }
   static void deleteArray_maplEintcOstringgR(void *p) {
      delete [] (static_cast<map<int,string>*>(p));
   }
   static void destruct_maplEintcOstringgR(void *p) {
      typedef map<int,string> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<int,string>

namespace ROOT {
   static TClass *maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR_Dictionary();
   static void maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR_TClassManip(TClass*);
   static void *new_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p = nullptr);
   static void *newArray_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(Long_t size, void *p);
   static void delete_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p);
   static void deleteArray_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p);
   static void destruct_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,shared_ptr<Pythia8::ParticleDataEntry> >*)
   {
      map<int,shared_ptr<Pythia8::ParticleDataEntry> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,shared_ptr<Pythia8::ParticleDataEntry> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,shared_ptr<Pythia8::ParticleDataEntry> >", -2, "map", 102,
                  typeid(map<int,shared_ptr<Pythia8::ParticleDataEntry> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,shared_ptr<Pythia8::ParticleDataEntry> >) );
      instance.SetNew(&new_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR);
      instance.SetNewArray(&newArray_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR);
      instance.SetDelete(&delete_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR);
      instance.SetDestructor(&destruct_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,shared_ptr<Pythia8::ParticleDataEntry> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<int,shared_ptr<Pythia8::ParticleDataEntry> >","std::map<int, std::shared_ptr<Pythia8::ParticleDataEntry>, std::less<int>, std::allocator<std::pair<int const, std::shared_ptr<Pythia8::ParticleDataEntry> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<int,shared_ptr<Pythia8::ParticleDataEntry> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<int,shared_ptr<Pythia8::ParticleDataEntry> >*>(nullptr))->GetClass();
      maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,shared_ptr<Pythia8::ParticleDataEntry> > : new map<int,shared_ptr<Pythia8::ParticleDataEntry> >;
   }
   static void *newArray_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,shared_ptr<Pythia8::ParticleDataEntry> >[nElements] : new map<int,shared_ptr<Pythia8::ParticleDataEntry> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p) {
      delete (static_cast<map<int,shared_ptr<Pythia8::ParticleDataEntry> >*>(p));
   }
   static void deleteArray_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p) {
      delete [] (static_cast<map<int,shared_ptr<Pythia8::ParticleDataEntry> >*>(p));
   }
   static void destruct_maplEintcOshared_ptrlEPythia8cLcLParticleDataEntrygRsPgR(void *p) {
      typedef map<int,shared_ptr<Pythia8::ParticleDataEntry> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<int,shared_ptr<Pythia8::ParticleDataEntry> >

namespace ROOT {
   static TClass *maplEintcOlonggR_Dictionary();
   static void maplEintcOlonggR_TClassManip(TClass*);
   static void *new_maplEintcOlonggR(void *p = nullptr);
   static void *newArray_maplEintcOlonggR(Long_t size, void *p);
   static void delete_maplEintcOlonggR(void *p);
   static void deleteArray_maplEintcOlonggR(void *p);
   static void destruct_maplEintcOlonggR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,long>*)
   {
      map<int,long> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,long>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,long>", -2, "map", 102,
                  typeid(map<int,long>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOlonggR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,long>) );
      instance.SetNew(&new_maplEintcOlonggR);
      instance.SetNewArray(&newArray_maplEintcOlonggR);
      instance.SetDelete(&delete_maplEintcOlonggR);
      instance.SetDeleteArray(&deleteArray_maplEintcOlonggR);
      instance.SetDestructor(&destruct_maplEintcOlonggR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,long> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<int,long>","std::map<int, long, std::less<int>, std::allocator<std::pair<int const, long> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<int,long>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOlonggR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<int,long>*>(nullptr))->GetClass();
      maplEintcOlonggR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOlonggR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOlonggR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,long> : new map<int,long>;
   }
   static void *newArray_maplEintcOlonggR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,long>[nElements] : new map<int,long>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOlonggR(void *p) {
      delete (static_cast<map<int,long>*>(p));
   }
   static void deleteArray_maplEintcOlonggR(void *p) {
      delete [] (static_cast<map<int,long>*>(p));
   }
   static void destruct_maplEintcOlonggR(void *p) {
      typedef map<int,long> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<int,long>

namespace ROOT {
   static TClass *maplEintcOdoublegR_Dictionary();
   static void maplEintcOdoublegR_TClassManip(TClass*);
   static void *new_maplEintcOdoublegR(void *p = nullptr);
   static void *newArray_maplEintcOdoublegR(Long_t size, void *p);
   static void delete_maplEintcOdoublegR(void *p);
   static void deleteArray_maplEintcOdoublegR(void *p);
   static void destruct_maplEintcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,double>*)
   {
      map<int,double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,double>", -2, "map", 102,
                  typeid(map<int,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,double>) );
      instance.SetNew(&new_maplEintcOdoublegR);
      instance.SetNewArray(&newArray_maplEintcOdoublegR);
      instance.SetDelete(&delete_maplEintcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEintcOdoublegR);
      instance.SetDestructor(&destruct_maplEintcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("map<int,double>","std::map<int, double, std::less<int>, std::allocator<std::pair<int const, double> > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const map<int,double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const map<int,double>*>(nullptr))->GetClass();
      maplEintcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,double> : new map<int,double>;
   }
   static void *newArray_maplEintcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) map<int,double>[nElements] : new map<int,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOdoublegR(void *p) {
      delete (static_cast<map<int,double>*>(p));
   }
   static void deleteArray_maplEintcOdoublegR(void *p) {
      delete [] (static_cast<map<int,double>*>(p));
   }
   static void destruct_maplEintcOdoublegR(void *p) {
      typedef map<int,double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class map<int,double>

namespace {
  void TriggerDictionaryInitialization_CreateDatasetDict_Impl() {
    static const char* headers[] = {
"CreateDataset.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/ktas/ge86rim/alice3-ca/sw/slc9_x86-64/pythia/v8311-28/include",
"/home/ktas/ge86rim/alice3-ca/sw/slc9_x86-64/ROOT/v6-32-06-alice1-17/include",
"/scratch5/ge86rim/alice3-ca/sw/slc9_x86-64/ROOT/v6-32-06-alice1-17/include/",
"/home/ktas/ge86rim/phsw/yaffa/yaffa/sim/pythia/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "CreateDatasetDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Basics.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Vec4;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Basics.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  RndmEngine;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Basics.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Rndm;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Info.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Info;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Settings.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Settings;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/StandardModel.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  AlphaStrong;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/StandardModel.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  AlphaEM;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/StandardModel.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  CoupSM;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ParticleData.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  DecayChannel;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ParticleData.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ParticleDataEntry;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ParticleData.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ParticleData;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Event.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Particle;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Event.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Junction;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Event.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Event;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceGeneric;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceGmZ;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceW;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceTop;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceFour;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceH;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceHchg;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceZprime;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceWprime;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceRhorizontal;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceExcited;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceGraviton;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceLeptoquark;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceNuRight;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceZRight;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceWRight;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceHchgchgLeft;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/ResonanceWidths.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  ResonanceHchgchgRight;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/SigmaProcess.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  InBeam;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/SigmaProcess.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  InPair;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/SigmaProcess.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  SigmaProcess;}
namespace Pythia8{class __attribute__((annotate("$clingAutoload$Pythia8/Pythia.h")))  __attribute__((annotate("$clingAutoload$CreateDataset.h")))  Pythia;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "CreateDatasetDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "CreateDataset.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Pythia8::AlphaEM", payloadCode, "@",
"Pythia8::AlphaStrong", payloadCode, "@",
"Pythia8::CoupSM", payloadCode, "@",
"Pythia8::DecayChannel", payloadCode, "@",
"Pythia8::Event", payloadCode, "@",
"Pythia8::InBeam", payloadCode, "@",
"Pythia8::InPair", payloadCode, "@",
"Pythia8::Info", payloadCode, "@",
"Pythia8::Junction", payloadCode, "@",
"Pythia8::Particle", payloadCode, "@",
"Pythia8::ParticleData", payloadCode, "@",
"Pythia8::ParticleDataEntry", payloadCode, "@",
"Pythia8::Pythia", payloadCode, "@",
"Pythia8::ResonanceExcited", payloadCode, "@",
"Pythia8::ResonanceFour", payloadCode, "@",
"Pythia8::ResonanceGeneric", payloadCode, "@",
"Pythia8::ResonanceGmZ", payloadCode, "@",
"Pythia8::ResonanceGraviton", payloadCode, "@",
"Pythia8::ResonanceH", payloadCode, "@",
"Pythia8::ResonanceHchg", payloadCode, "@",
"Pythia8::ResonanceHchgchgLeft", payloadCode, "@",
"Pythia8::ResonanceHchgchgRight", payloadCode, "@",
"Pythia8::ResonanceLeptoquark", payloadCode, "@",
"Pythia8::ResonanceNuRight", payloadCode, "@",
"Pythia8::ResonanceRhorizontal", payloadCode, "@",
"Pythia8::ResonanceTop", payloadCode, "@",
"Pythia8::ResonanceW", payloadCode, "@",
"Pythia8::ResonanceWRight", payloadCode, "@",
"Pythia8::ResonanceWprime", payloadCode, "@",
"Pythia8::ResonanceZRight", payloadCode, "@",
"Pythia8::ResonanceZprime", payloadCode, "@",
"Pythia8::Rndm", payloadCode, "@",
"Pythia8::RndmEngine", payloadCode, "@",
"Pythia8::Settings", payloadCode, "@",
"Pythia8::SigmaProcess", payloadCode, "@",
"Pythia8::Vec4", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("CreateDatasetDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_CreateDatasetDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_CreateDatasetDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_CreateDatasetDict() {
  TriggerDictionaryInitialization_CreateDatasetDict_Impl();
}
