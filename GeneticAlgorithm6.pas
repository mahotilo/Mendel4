{**************************************}
{*                                    *}
{*   ������ ������������� ���������   *}
{*           ������  4.7              *}
{*                                    *}
{**************************************}
{* ���������� ������                  *}
{* - ����������� ������������         *}
{* - ������� ���������, ����������    *}
{*   ���� ����������, ��������������� *}
{*   ������������� ����������� ���    *}
{*   ��������� ���������              *}
{* - ����������� ����������           *}
{*  (���� ��������� �������� +        *}
{*   ������ ��������� ��� �����)      *}
{* - ������ ����� ����������� ����    *}
{*   ���������                        *}
{**************************************}
{* ������������ ������                *}
{* - ����������� ������������         *}
{*   ������ �����������������         *}
{*   ������ � ���� �����              *}
{**************************************}
{* �����������                        *}
{* - �������� 62 ���� �� ����������   *}
{**************************************}
{* �� ���������                       *}
{* ��� ���������:                     *}
{*  ����������� �������� ������ - 0.1 *}
{*  �������� ��������� - 100%         *}
{* ��� �����:                         *}
{*  ������ ����������������� - �����  *}
{*  ����������� - ���                 *}
{*  ���������� - ���                  *}
{* ��� ���������:                     *}
{*  ����������� ������������          *}
{*  ����� ��������� - 0.5             *}
{*  ������� ������� - 1/250           *}
{*  ����������� ������������ - 0.005  *}
{*  ����������� �������� - 0.002      *}
{*  ����������� ����� ����������� - 1 *}
{*                                    *}
{* ����� ����������� � ������� - 2    *}
{**************************************}
{* 1998-2018 (c) ��������             *}
{*               ����������           *}
{*               ������������         *}
{*  Web: http://mahotilo.narod.ru     *}
{*  Email: mahotilo@narod.ru          *}
{**************************************}
{
31.01.18  v4.7g
 - ��� ����� ��������� � ���������� ����� �������� ����� ��������� 0
25.05.12  v4.7f
 * ��������� ���������� GAThreadsPerProcessor, ����������� ��������
  ���������� ������� �� ���������
19.02.11  v4.7e
 * � ������� RandomGamma ��������� ������ �� ������ �������� ln(0)
 * � ������� DaethProbability ��������� ������ �� ��������� �������� ln(0)
20.07.10  v4.7d
 * �������� �������������� ������������� ������ ������� �������
   � ������������� ������ �����
 * �������� BestIndividual �������� �� ���������� ��� �����������
   ����������� ���������� ��� ��������� � ��� � �������� ����������
25.05.10  v4.7�
 - ��������� ��� ������������� ������ ������� ������� � �������������
   ������ �����
15.11.09  v4.7b
 * �� TIndividual ������ ����� ����� ��� ���� ������ ���������� -
   ���������� �� Inf,Sup,Map � ����� �������� � ��������.
   ������ ��� ������ ��������� �� ������������ TGeneticAlgorithm.
   ��������� ������������� ������.
 * � TGeneticAlgorithm ������� ����� �������� � ������� ����� ���������
   � ����������, ����� �� ������������ �� � ������ ������ ���.
   +1% ��������������
 * � TGeneticAlgorithm ������� ����� ��� ���� ������ ��������� ���������
   +10% ��������������
 - ������ ����������� � TData_to_Bin, ����������� ��� ������������
   ������������ ��������� � ������� ������ ���� (> 53 ���)
11.11.09  v4.7�
 * �������� ����� �������� � ������������. ����������� ����������
   �������� ����������� ��� ������ ������ ���������
 * �������� �������� ����������� �������� � ������������ �� ���������
31.10.09  v4.7
 + ��� ���������� ������� ������� � ������������ ����� ������ �����������
   ��������/�������� ����� ������ �� ���������� ���
 + ����� ������� ��� ������������� ���������� ������� ������� ����� �����
   ����������� � ������� (�� ��������� 2)
 * �������� ����� ������������. ������ ������������ ������ �  ������ ����
   ���������� ���������� ����� ��������� ����� ����� ���������� ���������
28.08.09  v4.6a
 + ���������� ����������� � ��������� � ����� ������� - ��������� ��������
   � ������� ������������� �������� (������������� �����-�������������)
07.08.09  v4.6
 + ����� �������� ����� �� ���������� ������� ������ ���������� ������������
   � ������������ � ������ ���������
15.05.09  v4.5b
 + ������ �� �������� ����� ����� ��� ������� ���������� ������� �����
12.04.08  v4.5a
 + ��� ��������� �������� ������� ��������� inline (�������) ��� ���������
   �������� � Delphi 2007
04.04.08  v4.5
 + ���������� ������� ����� ������� �����
 * �������� ����� ������� ���� � ����������������� �� ������������ �������
   ���������.
 * ��������� ������������� ������������� ������� ��������.
 + ���������������� ����������� ����� �������� �� ��������� � �������
   TBitVector.CopyBitsTo()
09.12.07  v4.4e
 + ������ �� �������� Inf >= Sup ��� �������� ��
31.10.07  v4.4d
 * �������������� ������� �������� ���� ���� � ��������, ������� �������
   �� ����������
26.10.07  v4.4c
 + ����� ����������� ���������� ����� �������� �� ������ ��������, �� �����
   ����������� ������ � �� ���� �������� ������� ����������������� �����
30.08.07  v4.4b
 * ��������� ����������� �������������� ����������
22.07.07  v4.4a
 + ������������� ���������� ������� ����� �� �����, ����������� � ������ ������
07.04.07  v4.4
 + ������� ������� ������� �������� ������� � ������������ ������ ���
   ���������� ����������
07.06.06
09.04.06
12.01.06
20.02.05
28.03.02
19.11.99
02.11.99
22.04.99
31.01.99
20.08.98
25.02.98
08.02.98
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ����� ���������� �������������� ������ ���������� �������
// ����������������� ���������� ������ �������� ����� �����������
// � ������� (���������� ���������� GAProcessorsCount)
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

unit GeneticAlgorithm6;

interface

uses Classes,SysUtils,BaseMath,BitArray,SyncObjs;
const
 GAPopVer: byte = 47; //������ ������� ���������� ���������
 SupportedVersion: set of byte = [40,42,43,44,45,46];

type
  TMaxMin = (tcMinimize,tcMaximize);
  TEvalCondition = (ecNew,ecAll,ecNow);
  TOntogenesis = (ocNoOnt,ocOntAll);
  TDeath = (dcNoDeath,dcDeath);

  { Fitness function }
  TFitnFun = function(parameters: pointer; var condition: boolean): TData;

  TGeneticAlgorithm = class;

  TIndividual = class
//  private
   Pop: TGeneticAlgorithm;
   FChrLg,
   FChrRg,
   FChrLd,
   FChrRd: TBitVector;
   FAdge: integer;
   FParameters: TVector;
   FFitness: TData;
   FCondition: boolean;
   procedure MakeSexChromosome(ChrSg,ChrSd: TBitVector; CrMean,MtRate,TrProb,InProb: TData);
   procedure Crossover(ChrSg,ChrSd: TBitVector; CrMean: TData);  inline;
   procedure Translocation(ChrSg,ChrSd: TBitVector; TrProb: TData);  inline;
   procedure Inversion(ChrSg,ChrSd: TBitVector; InProb: TData);  inline;
   procedure Mutation(ChrSg,ChrSd: TBitVector; MtRate: TData);  inline;
   procedure DecodePhenotype;  inline;
   procedure Ontogenesis;
  public
   CrossoversCount,
   MutationsCount,
   InversionsCount,
   TranslocationsCount: integer;
   property ChrLg: TBitVector read FChrLg;
   property ChrLd: TBitVector read FChrLd;
   property ChrRg: TBitVector read FChrRg;
   property ChrRd: TBitVector read FChrRd;
   property Parameters: TVector read FParameters;
   property Fitness: TData read FFitness;
   property Condition: boolean read FCondition;
   property Adge: integer read FAdge;
   constructor Create(var GA: TGeneticAlgorithm);
   constructor LoadFromStream(Stream: TStream; var GA: TGeneticAlgorithm);
   procedure  SaveToStream(Stream: TStream);
   procedure SetRandom;
   procedure SetKnown(var Par: array of TData);
   destructor Destroy; override;
  end;


//!��� ����������� �������� ������������� ����� �������� ����� ��������� � ����� ������
  TStatistics = record
   CrossoversAv,
   MutationsAv,
   InversionsAv,
   TranslocationsAv: single;
   Epochs,
   Calculations,
   ParentsCount,
   ChildrensCount,
   DeadCount: integer;
   EfRadius,
   DifFitn: single;
   CreationTime,              // moment of population creation
   InitiationTime,            // calculation time for populationin itialization
   EvolutionTime: TDateTime;  // calculation time without time for population initialization
   NewBest: boolean;
   MeanAge: single;
   EpochTime: TDateTime;      // calculation time with correction for hibernation period
  end;


  TEvaluationThread = class(TThread)
  private
   FExecEvent: TEvent;
   Individual: TIndividual;
   FitnessFunction: TFitnFun;
  protected
   procedure Execute; override;
  public
   Runing: boolean;
   RunError: boolean;
   RunErrorMsg: string;
   procedure Run(Ind: TIndividual);
   constructor Create(FFunc: TFitnFun);
   destructor Destroy; override;
   procedure Terminate;
  end;


  TGeneticAlgorithm = class
  private
   FPopulation: TList;
   FPopulationSize: integer;
   FParametersVectorLen: integer;
   FChromosomeMap,
   FLinkGroupShift: TIntVector;
   FInfimum,
   FSupremum: TVector;
   FChromosomeLen: integer;
   FStatistics: TStatistics;
   FParentsCount,
   FChildrensCount,
   FDeadCount: integer;
   FStatCoeff: single;
   FitnessFunction: TFitnFun;
   FEvalThreadsPool: array of TEvaluationThread;
   CurVer: byte;
   TempCh,TempChg,TempChd: TBitVector;
   procedure AgeControl;
   procedure EvaluateAll;
   procedure MakeOffsprings;
   procedure EvaluateIndividual(Individual: TIndividual); inline;
   procedure EvaluateIndividualInThread(Individual: TIndividual); inline;
   procedure WaitForEvaluationThreads;
   procedure OntAll;
   procedure CalcStatistics;
   procedure ClearGenerationStatistics;
   function PopulationValid: boolean;
   function GetBestSolution: pointer;
   procedure SetUseEvaluationThreads(UET: boolean);
   procedure CalcChromosomeParams;
  public
   MaxMinTask: TMaxMin;
   ReevaluateCondition: TEvalCondition;
   OntogenesisCondition: TOntogenesis;
   DeathCondition: TDeath;
   MeanLifeLengh: integer;
   SelectionPressure,
   ActiveParents,
   CrossCountMean,
   MutationRate,
   TranslocationProbability,
   InversionProbability: TData;
   FUseEvaluationThreads: boolean;
   HibernationStartTime: TDateTime;
   HibernationStopTime: TDateTime;
   BestIndividual: TIndividual;
   property Population: TList read FPopulation;
   property PopulationSize: integer read FPopulationSize;
   property ParametersVectorLen: integer read FParametersVectorLen;
   property ChromosomeLen: integer read FChromosomeLen;
   property ChromosomeMap: TIntVector read FChromosomeMap;
   property Infimum: TVector read FInfimum;
   property Supremum: TVector read FSupremum;
   property Statistics: TStatistics read FStatistics;
   property ParentsCount: integer read FParentsCount;
   property ChildrensCount: integer read FChildrensCount;
   property BestSolution: pointer read GetBestSolution;
   property UseEvaluationThreads: boolean read FUseEvaluationThreads write SetUseEvaluationThreads;
   constructor Create(Dim: integer; var Inf,Sup,Map; Size: integer; MaxMin: TMaxMin);
   constructor Load(Stream: TStream);
   procedure Save(Stream: TStream);
   destructor Destroy; override;
   procedure SetFitnessFunction(FFunc: TFitnFun);
   procedure SetSelectionParametrs;
   procedure AddRandomToPopulation;
   procedure AddKnownToPopulation(var Par: array of TData);
   procedure DeleteFromPopulation(Index: integer);
   procedure Sort;
   procedure DoNextGeneration;
  end;


  procedure GetPredefGAParams(var SlP,APr,CrM,MtR,TrP,IvP: TData);
  function GetDefMeanLifeLengh(PopSize: integer): integer;
  procedure Dip_to_Gap(Dip1g,Dip1d,Dip2g,Dip2d,gap: TBitVector);
  procedure raiseError(const Msg: string); // !����� ������ ����������� �� inline ��� ���. ��������

var
 GAProcessorsCount,
 GAThreadsPerProcessor,
 GAThreadsCount: Integer;


implementation

var
 GAEvalThreadsCount: integer;
 FCSection: TCriticalSection;


procedure raiseError(const Msg: string); // ��� ��������� ��������� ���������� Exception ��������
begin
 raise Exception.Create(Msg);
end;


{-------------------------------- Math Tools ---------------------------------}
function RandomGamma(Mean: double): integer; inline;
{�����-�������������, ����� - ������������� ��������}
var
  tr: Extended;
  i,n: Integer;
begin
 tr := 1;
 n := trunc(Mean);
 if random <= Frac(Mean) then inc(n);
 for i := 1 to n do
 tr := tr * random;
 tr := -ln(tr + 1e-1000);
 if n = 1 then tr := tr*1.05;
 Result := round(tr);
end;


function Distance(v1,v2: TVector): TData; inline;
var i,h: integer;
begin
 result := 0;
 i := 0;
 h := high(v1);
 repeat
  result := result + sqr(v1[i]-v2[i]);
  inc(i);
 until i > h;
 result := sqrt(result);
end;


{-------------------------------- Bin Tools ---------------------------------}
procedure Bin_to_Gray(bin, gray: TBitVector; index,count: integer);  inline;
var
 i: integer;
 x,y: byte;
begin
 x := 0;
 for i := index+count-1 downto index do
 begin
  if bin[i] then y := 1 else y := 0;
  gray[i] := (x xor y) = 1;
  x := y;
 end;
end;


procedure Gray_to_Bin(gray, bin : TBitVector; index,count: integer); inline;
var
 i: integer;
 x: byte;
begin
 x := 0;
 for i := index+count-1 downto index do
 begin
  if gray[i] then x := 1 - x;
  bin[i] := x = 1;
 end;
end;


procedure Bin_to_TData(bin: TBitVector; min,max: TData; var Data: TData; index,count: integer); inline;
var
 d: Int64;
 i: integer;
begin
 Data := 0;
 d := 1;
 for i := index to index+count-1 do
 begin
  if bin[i] then Data := Data + d;
  d := d shl 1;
 end;
 Data := Data*(max-min)/(d-1) + min;
end;


procedure TData_to_Bin(rl,min,max: TData; bin: TBitVector;  index,count: integer); inline;
var
 rle: Extended;
 t,rli: Int64;
 i: integer;
begin
 t := 1;
 t := t shl count;
 rle := (rl-min)/(max-min)*(t-1);
 rli := trunc(rle+0.5);
 for i := index+count-1 downto index  do
 begin
  t := t shr 1;
  bin[i] := (rli and t) > 0;
 end;
end;


{-------------------------------- Gen Tools ---------------------------------}
procedure Dip_to_Gap(Dip1g,Dip1d,Dip2g,Dip2d,gap: TBitVector); inline;
var
 i: integer;
begin
 i := 1;
 repeat // for i := 1 to gap.Size do
  if Dip1g[i] and Dip1d[i] then gap[i] := true else       {Dx}
  if not Dip1g[i] and Dip1d[i] then gap[i] := false else  {dx}
  if Dip2g[i]                                             {xD or xR}
  then gap[i] := true
  else gap[i] := false;
  inc(i);
 until i > gap.Size;
end;



{-------------------------------- TIndividual ---------------------------------}
constructor TIndividual.Create(var GA: TGeneticAlgorithm);
begin
 Pop := GA;
 inherited Create;

 FAdge := 0;
 FChrLg := TBitVector.Create;
 FChrLd := TBitVector.Create;
 FChrRg := TBitVector.Create;
 FChrRd := TBitVector.Create;
 FChrLg.Size := Pop.FChromosomeLen;
 FChrLd.Size := Pop.FChromosomeLen;
 FChrRg.Size := Pop.FChromosomeLen;
 FChrRd.Size := Pop.FChromosomeLen;
 SetLength(FParameters,Pop.FParametersVectorLen);
end;


constructor TIndividual.LoadFromStream(Stream: TStream; var GA: TGeneticAlgorithm);
var Dummy: integer;
begin
 inherited Create;
 Pop := GA;

{�������������}
 if Pop.CurVer < 47 then
 begin
  Stream.ReadBuffer(Dummy, SizeOf(Integer));
  Stream.ReadBuffer(Dummy, SizeOf(Integer));
 end;
{^^^^^^^^^^^^^}
 Stream.ReadBuffer(FFitness, SizeOf(FFitness));
 Stream.ReadBuffer(FCondition, SizeOf(FCondition));
 Stream.ReadBuffer(CrossoversCount, SizeOf(CrossoversCount));
 Stream.ReadBuffer(MutationsCount, SizeOf(MutationsCount));
 Stream.ReadBuffer(InversionsCount, SizeOf(InversionsCount));
 Stream.ReadBuffer(TranslocationsCount, SizeOf(TranslocationsCount));
 Stream.ReadBuffer(FAdge, SizeOf(FAdge));
 FChrLg := TBitVector.LoadFromStream(Stream);
 FChrLd := TBitVector.LoadFromStream(Stream);
 FChrRg := TBitVector.LoadFromStream(Stream);
 FChrRd := TBitVector.LoadFromStream(Stream);
 LoadTVectorFromStream(Stream,FParameters);
end;


procedure TIndividual.SaveToStream(Stream: TStream);
begin
 Stream.WriteBuffer(FFitness, SizeOf(FFitness));
 Stream.WriteBuffer(FCondition, SizeOf(FCondition));
 Stream.WriteBuffer(CrossoversCount, SizeOf(CrossoversCount));
 Stream.WriteBuffer(MutationsCount, SizeOf(MutationsCount));
 Stream.WriteBuffer(InversionsCount, SizeOf(InversionsCount));
 Stream.WriteBuffer(TranslocationsCount, SizeOf(TranslocationsCount));
 Stream.WriteBuffer(FAdge, SizeOf(FAdge));
 FChrLg.SaveToStream(Stream);
 FChrLd.SaveToStream(Stream);
 FChrRg.SaveToStream(Stream);
 FChrRd.SaveToStream(Stream);
 SaveTVectorToStream(Stream,FParameters);
end;


destructor TIndividual.Destroy;
begin
 Finalize(FParameters);
 FChrLg.Free;
 FChrLd.Free;
 FChrRg.Free;
 FChrRd.Free;
 inherited Destroy;
end;


procedure TIndividual.SetRandom;
var i,r: integer;
begin
 for i := 1 to Pop.FChromosomeLen do
 begin
  r := random(4)+1;
  case r of
  1: begin FChrLg[i] := true; FChrLd[i] := true; end;   {D}
  2: begin FChrLg[i] := true; FChrLd[i] := false; end;  {R}
  3: begin FChrLg[i] := false; FChrLd[i] := false; end; {r}
  4: begin FChrLg[i] := false; FChrLd[i] := true; end;  {d}
  end;

  r := random(4)+1;
  case r of
  1: begin FChrRg[i] := true; FChrRd[i] := true; end;   {D}
  2: begin FChrRg[i] := true; FChrRd[i] := false; end;  {R}
  3: begin FChrRg[i] := false; FChrRd[i] := false; end; {r}
  4: begin FChrRg[i] := false; FChrRd[i] := true; end;  {d}
  end;
 end;
 DecodePhenotype;
end;


procedure TIndividual.SetKnown(var Par: array of TData);
var
 i: integer;
begin
 move(Par[0],FParameters[0],Sizeof(TData)*Length(FParameters));
 for i := 0 to Pop.FParametersVectorLen-1 do
 begin
  TData_to_Bin(FParameters[i],Pop.FInfimum[i],Pop.FSupremum[i],Pop.TempCh,Pop.FLinkGroupShift[i]+1,Pop.FChromosomeMap[i]);
  Bin_to_Gray(Pop.TempCh,Pop.TempCh,Pop.FLinkGroupShift[i]+1,Pop.FChromosomeMap[i]);
 end;

 for i := 1 to Pop.FChromosomeLen do
 begin
  if Pop.TempCh[i]
  then
   begin
    if random < 0.5
    then
     begin
      FChrLg[i] := true; FChrLd[i] := true;   {D}
      FChrRg[i] := true; FChrRd[i] := true;   {D}
     end
    else
     begin
      FChrLg[i] := true; FChrLd[i] := false;  {R}
      FChrRg[i] := true; FChrRd[i] := false;  {R}
     end
   end
  else
   begin
    if random < 0.5
    then
     begin
      FChrLg[i] := false; FChrLd[i] := true;  {d}
      FChrRg[i] := false; FChrRd[i] := true;  {d}
     end
    else
     begin
      FChrLg[i] := false; FChrLd[i] := false; {r}
      FChrRg[i] := false; FChrRd[i] := false  {r}
     end
   end;
 end;

 DecodePhenotype;
end;


procedure TIndividual.Crossover(ChrSg,ChrSd: TBitVector; CrMean: TData);
var
 CrossPoint,ChrBeg,i,j: integer;
begin
 ChrSg.CopyBitVectorFrom(FChrLg);
 ChrSd.CopyBitVectorFrom(FChrLd);
 Pop.TempChg.CopyBitVectorFrom(FChrRg);
 Pop.TempChd.CopyBitVectorFrom(FChrRd);

 CrossPoint := 1;
 for i := 0 to Pop.FParametersVectorLen-1 do
 begin
  if random < 0.5 then
  begin
   FChrRg.CopyBitsTo(ChrSg,CrossPoint,Pop.FChromosomeMap[i],CrossPoint);
   FChrRd.CopyBitsTo(ChrSd,CrossPoint,Pop.FChromosomeMap[i],CrossPoint);
   FChrLg.CopyBitsTo(Pop.TempChg,CrossPoint,Pop.FChromosomeMap[i],CrossPoint);
   FChrLd.CopyBitsTo(Pop.TempChd,CrossPoint,Pop.FChromosomeMap[i],CrossPoint);
  end;
  CrossPoint := CrossPoint + Pop.FChromosomeMap[i];
 end;

 ChrBeg := 0;
 for i := 0 to Pop.FParametersVectorLen-1 do
 begin
  for j := 1 to RandomGamma(CrMean) do
  begin
   CrossPoint := random(Pop.FChromosomeMap[i]-1)+1;

   ChrSg.CopyBitsTo(Pop.TempCh,ChrBeg+CrossPoint+1,Pop.FChromosomeMap[i]-CrossPoint,ChrBeg+CrossPoint+1);
   Pop.TempChg.CopyBitsTo(ChrSg,ChrBeg+CrossPoint+1,Pop.FChromosomeMap[i]-CrossPoint,ChrBeg+CrossPoint+1);
   Pop.TempCh.CopyBitsTo(Pop.TempChg,ChrBeg+CrossPoint+1,Pop.FChromosomeMap[i]-CrossPoint,ChrBeg+CrossPoint+1);

   ChrSd.CopyBitsTo(Pop.TempCh,ChrBeg+CrossPoint+1,Pop.FChromosomeMap[i]-CrossPoint,ChrBeg+CrossPoint+1);
   Pop.TempChd.CopyBitsTo(ChrSd,ChrBeg+CrossPoint+1,Pop.FChromosomeMap[i]-CrossPoint,ChrBeg+CrossPoint+1);
   Pop.TempCh.CopyBitsTo(Pop.TempChd,ChrBeg+CrossPoint+1,Pop.FChromosomeMap[i]-CrossPoint,ChrBeg+CrossPoint+1);
   Inc(CrossoversCount);
  end;
  ChrBeg := ChrBeg + Pop.FChromosomeMap[i];
 end;

 if random > 0.5 then
 begin
  ChrSg.CopyBitVectorFrom(Pop.TempChg);
  ChrSd.CopyBitVectorFrom(Pop.TempChd);
 end;
end;


procedure TIndividual.Translocation(ChrSg,ChrSd: TBitVector; TrProb: TData);
var
 n,i,chn1,chn2,ChLen1,ChLen2,Len,TransPoint1,TransPoint2: integer;
begin
 if (Pop.FParametersVectorLen > 1) then
 for n := 1 to RandomGamma(Pop.FParametersVectorLen*TrProb) do
 begin
  chn1 := random(Pop.FParametersVectorLen);
  chn2 := random(Pop.FParametersVectorLen);
  if chn2 = chn1 then inc(chn2);
  if chn2 >= Pop.FParametersVectorLen then chn2 := 0;

  if Pop.FChromosomeMap[chn1] > Pop.FChromosomeMap[chn2] then
  begin
   i := chn1;
   chn1 := chn2;
   chn2 := i;
  end;

  ChLen1 := Pop.FChromosomeMap[chn1];
  ChLen2 := Pop.FChromosomeMap[chn2];
  if (ChLen1 < 2) or (ChLen2 < 2) then exit;

  Len := 1+random(ChLen1-1);
  if random < 0.5
  then
   begin // ������ ������
    TransPoint1 := 1;
    TransPoint2 := 1;
   end
  else
   begin // ������ ������
    TransPoint1 := ChLen1-Len+1;;
    TransPoint2 := ChLen2-Len+1;
   end;

  TransPoint1 := Pop.FLinkGroupShift[chn1] + TransPoint1;
  TransPoint2 := Pop.FLinkGroupShift[chn2] + TransPoint2;

  Pop.TempCh.CopyBitVectorFrom(ChrSg);
  Pop.TempCh.CopyBitsTo(ChrSg,TransPoint1,Len,TransPoint2);
  Pop.TempCh.CopyBitsTo(ChrSg,TransPoint2,Len,TransPoint1);

  Pop.TempCh.CopyBitVectorFrom(ChrSd);
  Pop.TempCh.CopyBitsTo(ChrSd,TransPoint1,Len,TransPoint2);
  Pop.TempCh.CopyBitsTo(ChrSd,TransPoint2,Len,TransPoint1);

  Inc(TranslocationsCount);
 end;
end;


procedure TIndividual.Inversion(ChrSg,ChrSd: TBitVector; InProb: TData);
var
 ChLen,InvPoint,Len,n,i,chn: integer;
 bit: boolean;
begin
 for n := 1 to RandomGamma(Pop.FParametersVectorLen*InProb) do
 begin
  chn := random(Pop.FParametersVectorLen);
  if (Pop.FChromosomeMap[chn] > 2) then
  begin
   ChLen := Pop.FChromosomeMap[chn];
   InvPoint := 1+random(ChLen-2);
   Len := random(ChLen-InvPoint-1) + 2 ;
   InvPoint := Pop.FLinkGroupShift[chn] + InvPoint;
   for i := 1 to Len div 2 do
   begin
    bit := ChrSg[InvPoint+i];
    ChrSg[InvPoint+i] := ChrSg[InvPoint+1+Len-i];
    ChrSg[InvPoint+1+Len-i] := bit;
    bit := ChrSd[InvPoint+i];
    ChrSd[InvPoint+i] := ChrSd[InvPoint+1+Len-i];
    ChrSd[InvPoint+1+Len-i] := bit;
   end;
   Inc(InversionsCount);
  end;
 end;
end;


procedure TIndividual.Mutation(ChrSg,ChrSd: TBitVector; MtRate: TData);
var
 MutPoint,i,n: integer;
begin
 for i := 1 to RandomGamma(Pop.FChromosomeLen*MtRate) do
 begin
  MutPoint := random(Pop.FChromosomeLen)+1;
  n := random(3); // 0 - g;  1 - d;  2 - gd
  if n <> 1 then
  ChrSg[MutPoint] := not ChrSg[MutPoint];
  if n <> 0 then
  ChrSd[MutPoint] := not ChrSd[MutPoint];
  Inc(MutationsCount);
 end;
end;


procedure TIndividual.MakeSexChromosome(ChrSg,ChrSd: TBitVector; CrMean,MtRate,TrProb,InProb: TData);
begin
 CrossoversCount := 0;
 TranslocationsCount := 0;
 InversionsCount := 0;
 MutationsCount := 0;
 Crossover(ChrSg,ChrSd,CrMean);
 Translocation(ChrSg,ChrSd,TrProb);
 Inversion(ChrSg,ChrSd,InProb);
 Mutation(ChrSg,ChrSd,MtRate);
end;


procedure TIndividual.Ontogenesis;
var
 i: integer;
begin
 for i := 1 to Pop.FChromosomeLen do
 if  (FChrLd[i] and FChrRd[i]) and (FChrLg[i] xor FChrRg[i]) then
 if random < 0.5 then
 if random < 0.5
 then FChrLd[i] := not FChrLd[i]
 else FChrRd[i] := not FChrRd[i];
 DecodePhenotype;
end;


procedure TIndividual.DecodePhenotype;
var
 i: integer;
begin
 Dip_to_Gap(FChrLg,FChrLd,FChrRg,FChrRd,Pop.TempCh);
 for i := 0 to Pop.FParametersVectorLen-1 do
 begin
  Gray_to_Bin(Pop.TempCh,Pop.TempCh,Pop.FLinkGroupShift[i]+1,Pop.FChromosomeMap[i]);
  Bin_to_TData(Pop.TempCh,Pop.FInfimum[i],Pop.FSupremum[i],FParameters[i],Pop.FLinkGroupShift[i]+1,Pop.FChromosomeMap[i]);
 end;
end;



{---------------------------- TEvaluationThread -------------------------------}
constructor TEvaluationThread.Create(FFunc: TFitnFun);
begin
 Runing := false;
 FExecEvent := TEvent.Create(false);
 FExecEvent.ResetEvent;
 FreeOnTerminate := True;
 FitnessFunction := FFunc;
 inherited Create(false);
end;


destructor TEvaluationThread.Destroy;
begin
 FExecEvent.Free;
 inherited Destroy;
end;


procedure TEvaluationThread.Terminate;
begin
 inherited Terminate;
 FExecEvent.SetEvent;
end;


procedure TEvaluationThread.Run(Ind: TIndividual);
begin
 Individual := Ind;
 RunError := false;
 FCSection.Acquire;
  Inc(GAEvalThreadsCount);
 FCSection.Release;
 Runing := true;
 FExecEvent.SetEvent;
end;


procedure TEvaluationThread.Execute;
begin
 repeat
  FExecEvent.WaitFor(INFINITE);
  if Runing then
  begin
   with Individual do
   try
    FCondition := true;
    FFitness := FitnessFunction(FParameters,FCondition);
   except
    on E: Exception do
    begin
     RunError := true;
     RunErrorMsg := E.Message;
    end
   end;
   FCSection.Acquire;
   Dec(GAEvalThreadsCount);
   FCSection.Release;
   FExecEvent.ResetEvent;
   Runing := false;
  end;
  sleep(0);
 until Terminated;
end;



{---------------------------- TGeneticAlgorithm -------------------------------}
procedure GetPredefGAParams(var SlP,APr,CrM,MtR,TrP,IvP: TData); inline;
begin
 SlP := 0.1;
 APr := 1;
 CrM := 1;
 MtR := 1/250;
 TrP := 0.005;
 IvP := 0.002;
end;


function GetDefMeanLifeLengh(PopSize: integer): integer; inline;
begin
 result := PopSize*100; {���� ��������� �������� - size*30 ��������� ��� ����������� ��������� ���� + ����������� ���� ������� ������-��������� ��������� ����}
 if result < 300 then result := 300;
end;



constructor TGeneticAlgorithm.Create(Dim: integer; var Inf,Sup,Map; Size: integer; MaxMin: TMaxMin);
var
 i: integer;
 FloatMap: TVector;
begin
 inherited Create;
 FStatistics.CreationTime := Now;
 FStatistics.EvolutionTime := 0;
 FStatistics.InitiationTime := 0;

 CurVer := GAPopVer;

 FitnessFunction := nil;
 FParametersVectorLen := Dim;

 SetLength(FInfimum,FParametersVectorLen);
 SetLength(FSupremum,FParametersVectorLen);
 move(Inf,FInfimum[0],sizeof(TData)*FParametersVectorLen);
 move(Sup,FSupremum[0],sizeof(TData)*FParametersVectorLen);

 for i := 0 to FParametersVectorLen-1 do
 if FInfimum[i] >= FSupremum[i] then
 raiseError('Infinum greater or equel to supremum');

 SetLength(FChromosomeMap,FParametersVectorLen);
 SetLength(FLinkGroupShift,FParametersVectorLen);
 SetLength(FloatMap,FParametersVectorLen);
 move(Map,FloatMap[0],sizeof(TData)*FParametersVectorLen);
 for i := 0 to high(FChromosomeMap) do
 FChromosomeMap[i] := round(FloatMap[i]);
 Finalize(FloatMap);
 CalcChromosomeParams;

 FStatistics.Epochs := 0;
 FStatistics.Calculations := 0;
 FStatistics.EpochTime := 0;
 FStatistics.NewBest := true;
 FPopulationSize := 0;
 MaxMinTask := MaxMin;
 ReevaluateCondition := ecNew;
 OntogenesisCondition := ocNoOnt;
 DeathCondition := dcNoDeath;
 MeanLifeLengh := 0;
 FDeadCount := 0;
 FUseEvaluationThreads := false;
 GAEvalThreadsCount := 0;
 HibernationStartTime := 0;
 HibernationStopTime := 0;

 GetPredefGAParams(SelectionPressure,ActiveParents,CrossCountMean,
                 MutationRate,TranslocationProbability,InversionProbability);

 TempCh := TBitVector.Create;
 TempCh.Size := FChromosomeLen;
 TempChg := TBitVector.Create;
 TempChg.Size := FChromosomeLen;
 TempChd := TBitVector.Create;
 TempChd.Size := FChromosomeLen;

 if Size < 2 then raiseError('Incorrect population size');
 try
  FPopulation := TList.Create;
  FPopulation.Capacity := Size;
 except
  on E: Exception do
   begin
    FPopulation := nil;
    E.Message := 'Out of memory.'+#13+'Can not create population';
    raise;
   end;
 end;
 BestIndividual := nil;

 FStatistics.InitiationTime := Now - FStatistics.CreationTime;
end;


constructor TGeneticAlgorithm.Load(Stream: TStream);
var
 i: integer;
 Ind: TIndividual;
 vererr,UET: boolean;
 Dummy: TData;
 Vect: TVector;
begin
 inherited Create;
 vererr := false;
 with Stream do
 try
  FitnessFunction := nil;
  GAEvalThreadsCount := 0;
  ReadBuffer(CurVer, SizeOf(CurVer));
  if CurVer <> GAPopVer then
  vererr := not (CurVer in SupportedVersion);
  if vererr then raiseError('Unsupported GA version');

  ReadBuffer(FParametersVectorLen, SizeOf(FParametersVectorLen));
  LoadTVectorFromStream(Stream,FInfimum);
  LoadTVectorFromStream(Stream,FSupremum);
{�������������}
  if CurVer < 47 then
  begin
   LoadTVectorFromStream(Stream,Vect);
   SetLength(FChromosomeMap,FParametersVectorLen);
   for i := 0 to high(Vect) do
   FChromosomeMap[i] := round(Vect[i]);
   Finalize(Vect);
  end
  else
{^^^^^^^^^^^^^}
  LoadTIntVectorFromStream(Stream,FChromosomeMap);
  SetLength(FLinkGroupShift,FParametersVectorLen);
  CalcChromosomeParams;
  ReadBuffer(FParentsCount, SizeOf(FParentsCount));
  ReadBuffer(FChildrensCount, SizeOf(FChildrensCount));
  ReadBuffer(FDeadCount, SizeOf(FDeadCount));
  ReadBuffer(SelectionPressure, SizeOf(SelectionPressure));
  ReadBuffer(ActiveParents, SizeOf(ActiveParents));
  ReadBuffer(CrossCountMean, SizeOf(CrossCountMean));
{�������������}
  if CurVer < 46 then
  ReadBuffer(Dummy, SizeOf(TData));
{^^^^^^^^^^^^^}
  ReadBuffer(MutationRate, SizeOf(MutationRate));
  ReadBuffer(TranslocationProbability, SizeOf(TranslocationProbability));
  ReadBuffer(InversionProbability, SizeOf(InversionProbability));
  ReadBuffer(FPopulationSize, SizeOf(FPopulationSize));
  ReadBuffer(MaxMinTask, SizeOf(MaxMinTask));
  ReadBuffer(ReevaluateCondition, SizeOf(ReevaluateCondition));
  ReadBuffer(OntogenesisCondition, SizeOf(OntogenesisCondition));
  ReadBuffer(DeathCondition, SizeOf(DeathCondition));
{�������������}
  if CurVer < 47 then
  GetPredefGAParams(SelectionPressure,ActiveParents,CrossCountMean,
                    MutationRate,TranslocationProbability,InversionProbability);
{^^^^^^^^^^^^^}
{�������������}
  if CurVer >= 44
  then
   begin
    ReadBuffer(UET, SizeOf(UseEvaluationThreads));
    UseEvaluationThreads := UET;
   end
  else
   FUseEvaluationThreads := false;
{^^^^^^^^^^^^^}
{�������������}
  case CurVer of
  40..44: ReadBuffer(FStatistics, 80);
  else
   ReadBuffer(FStatistics, SizeOf(FStatistics));
  end;
  if CurVer < 43 then
  begin
   FStatistics.InitiationTime := 0;
   FStatistics.EpochTime := 0;
  end;
{^^^^^^^^^^^^^}
  ReadBuffer(MeanLifeLengh, SizeOf(MeanLifeLengh));

  try
   FPopulation := TList.Create;
   FPopulation.Capacity := FPopulationSize;
  except
   on E: Exception do
    begin
     E.Message := 'Can not create population';
     raise;
    end;
  end;
  for i := 0 to FPopulationSize-1 do
  begin
   Ind := TIndividual.LoadFromStream(Stream,Self);
   FPopulation.Add(Ind);
  end;

  BestIndividual := TIndividual(Population[0]);
  TempCh := TBitVector.Create;
  TempCh.Size := FChromosomeLen;
  TempChg := TBitVector.Create;
  TempChg.Size := FChromosomeLen;
  TempChd := TBitVector.Create;
  TempChd.Size := FChromosomeLen;
 except
  on E: Exception do
  begin
   if vererr
   then E.Message := 'Incompatible version of population'
   else E.Message := 'Population load error';
   raise;
  end;
 end;
end;


procedure TGeneticAlgorithm.Save(Stream: TStream);
var i: integer;
begin
 with Stream do
 try
  WriteBuffer(GAPopVer, SizeOf(GAPopVer));
  WriteBuffer(FParametersVectorLen, SizeOf(FParametersVectorLen));
  SaveTVectorToStream(Stream,FInfimum);
  SaveTVectorToStream(Stream,FSupremum);
  SaveTIntVectorToStream(Stream,FChromosomeMap);
  WriteBuffer(FParentsCount, SizeOf(FParentsCount));
  WriteBuffer(FChildrensCount, SizeOf(FChildrensCount));
  WriteBuffer(FDeadCount, SizeOf(FDeadCount));
  WriteBuffer(SelectionPressure, SizeOf(SelectionPressure));
  WriteBuffer(ActiveParents, SizeOf(ActiveParents));
  WriteBuffer(CrossCountMean, SizeOf(CrossCountMean));
  WriteBuffer(MutationRate, SizeOf(MutationRate));
  WriteBuffer(TranslocationProbability, SizeOf(TranslocationProbability));
  WriteBuffer(InversionProbability, SizeOf(InversionProbability));
  WriteBuffer(FPopulationSize, SizeOf(FPopulationSize));
  WriteBuffer(MaxMinTask, SizeOf(MaxMinTask));
  WriteBuffer(ReevaluateCondition, SizeOf(ReevaluateCondition));
  WriteBuffer(OntogenesisCondition, SizeOf(OntogenesisCondition));
  WriteBuffer(DeathCondition, SizeOf(DeathCondition));
  WriteBuffer(UseEvaluationThreads, SizeOf(UseEvaluationThreads));
  WriteBuffer(FStatistics, SizeOf(FStatistics));
  WriteBuffer(MeanLifeLengh, SizeOf(MeanLifeLengh));
  for i := 0 to FPopulationSize-1 do
  TIndividual(FPopulation[i]).SaveToStream(Stream);
 except
  on E: Exception do
  begin
   E.Message := 'Population save error';
   raise;
  end;
 end;
end;


destructor TGeneticAlgorithm.Destroy;
var i: integer;
begin
 Finalize(FInfimum);
 Finalize(FSupremum);
 Finalize(FChromosomeMap);
 Finalize(FLinkGroupShift);
 if Assigned(FPopulation) then
 begin
  for i := 0 to FPopulation.Count-1 do
   TIndividual(FPopulation[i]).Free;
  FPopulation.Free;
 end;
 UseEvaluationThreads := false;
 TempCh.Free;
 TempChg.Free;
 TempChd.Free;
 inherited Destroy;
end;


procedure TGeneticAlgorithm.SetFitnessFunction(FFunc: TFitnFun);
var i: integer;
begin
 FitnessFunction := FFunc;
 if UseEvaluationThreads then
 for i := 0 to high(FEvalThreadsPool) do
 FEvalThreadsPool[i].FitnessFunction := FFunc;
end;


procedure TGeneticAlgorithm.SetSelectionParametrs;
begin
 FChildrensCount := round(FPopulationSize*SelectionPressure);
 if DeathCondition = dcDeath then
 begin
  if FChildrensCount < FDeadCount then FChildrensCount := FDeadCount;
  if FChildrensCount > FPopulationSize - 2 then FChildrensCount := FPopulationSize - 2;
 end;
 if FChildrensCount <= 0 then FChildrensCount := 1;
 FParentsCount := round((FPopulationSize - FChildrensCount)*ActiveParents);
 if FParentsCount = 0 then FParentsCount := 2;
 FStatCoeff := 0.5/FChildrensCount;
end;


function TGeneticAlgorithm.PopulationValid: boolean;
begin
 result := true;
 result := result and (FPopulationSize > 1);
 result := result and Assigned(FitnessFunction);
end;


procedure TGeneticAlgorithm.SetUseEvaluationThreads(UET: boolean);
var i: integer;
begin
 if not FUseEvaluationThreads and UET then
 begin
  GAThreadsCount := GAProcessorsCount*GAThreadsPerProcessor;
  SetLength(FEvalThreadsPool,GAThreadsCount);
  for i := 0 to high(FEvalThreadsPool) do
  FEvalThreadsPool[i] := TEvaluationThread.Create(FitnessFunction);
 end;
 if FUseEvaluationThreads and not UET then
 begin
  for i := 0 to high(FEvalThreadsPool) do
  FEvalThreadsPool[i].Terminate;
  SetLength(FEvalThreadsPool,0);
 end;
 FUseEvaluationThreads := UET;
end;


procedure TGeneticAlgorithm.AddRandomToPopulation;
var
 Ind: TIndividual;
 StartTime: TDateTime;
begin
 if not Assigned(FitnessFunction) then raiseError('No fitness function defined');
 StartTime := Now;
 Ind := TIndividual.Create(Self);
 Ind.SetRandom;
 EvaluateIndividual(Ind);
 FPopulation.Add(Ind);
 Inc(FPopulationSize);
 Inc(FStatistics.Calculations);
 if FPopulationSize > 1 then
 begin
  SetSelectionParametrs;
  MeanLifeLengh := GetDefMeanLifeLengh(FPopulationSize);
  Sort;
  CalcStatistics;
 end;
 with FStatistics do InitiationTime := InitiationTime + Now - StartTime;
end;


procedure TGeneticAlgorithm.AddKnownToPopulation(var Par: array of TData);
var
 Ind: TIndividual;
 StartTime: TDateTime;
begin
 if not Assigned(FitnessFunction) then raiseError('No fitness function defined');
 StartTime := Now;
 Ind := TIndividual.Create(Self);
 Ind.SetKnown(Par);
 EvaluateIndividual(Ind);
 FPopulation.Add(Ind);
 Inc(FPopulationSize);
 Inc(FStatistics.Calculations);
 if FPopulationSize > 1 then
 begin
  SetSelectionParametrs;
  MeanLifeLengh := GetDefMeanLifeLengh(FPopulationSize);
  Sort;
  CalcStatistics;
 end;
 with FStatistics do InitiationTime := InitiationTime + Now - StartTime;
end;


procedure TGeneticAlgorithm.DeleteFromPopulation(Index: integer);
begin
 if FPopulation.Count > 0  then
 begin
  TIndividual(FPopulation[Index]).Free;
  FPopulation.Delete(Index);
  FPopulation.Capacity := FPopulation.Count;
  Dec(FPopulationSize);
  SetSelectionParametrs;
  MeanLifeLengh := GetDefMeanLifeLengh(FPopulationSize);
  CalcStatistics;
 end;
end;


procedure TGeneticAlgorithm.EvaluateIndividual(Individual: TIndividual);
begin
 with Individual do
 begin
  FCondition := true;
  FFitness := FitnessFunction(FParameters,FCondition);
 end;
end;


procedure TGeneticAlgorithm.EvaluateIndividualInThread(Individual: TIndividual);
var
 i: integer;
begin
 repeat
  i := 0;
  while (i <= high(FEvalThreadsPool)) and FEvalThreadsPool[i].Runing do
  inc(i);
  sleep(0);
 until (i <= high(FEvalThreadsPool));
 if not FEvalThreadsPool[i].RunError then
 FEvalThreadsPool[i].Run(Individual);
end;


procedure TGeneticAlgorithm.WaitForEvaluationThreads;
var
 i: integer;
begin
 while GAEvalThreadsCount > 0 do
 sleep(0);
 for i := 0 to high(FEvalThreadsPool) do
 if FEvalThreadsPool[i].RunError then
 raiseError(FEvalThreadsPool[i].RunErrorMsg);
end;


function TGeneticAlgorithm.GetBestSolution: pointer;
begin
 result := BestIndividual.Parameters;
end;


procedure TGeneticAlgorithm.EvaluateAll;
var i,n: integer;
begin
 if ReevaluateCondition = ecAll
 then n := FParentsCount
 else n := FPopulationSize; // ecNow
 try
  for i := 0 to n-1 do
  begin
   if UseEvaluationThreads
   then EvaluateIndividualInThread(TIndividual(FPopulation[i]))
   else EvaluateIndividual(TIndividual(FPopulation[i]));
   Inc(FStatistics.Calculations);
  end;
 finally
  if UseEvaluationThreads
  then WaitForEvaluationThreads;
 end;
 FStatistics.NewBest := true;
end;


procedure TGeneticAlgorithm.OntAll;
var i: integer;
    Ind: TIndividual;
    OldFit: TData;
    OldP: TVector;
begin
 for i := 1 to FParentsCount do
 if random < 0.3 then
 begin
  Ind := TIndividual(FPopulation[i-1]);
  OldFit := Ind.Fitness;
  OldP := Copy(Ind.Parameters);
  Ind.Ontogenesis;
  EvaluateIndividual(Ind);
  if (Ind.Fitness > OldFit) xor (MaxMinTask = tcMaximize)
  then
   begin
    Ind.FFitness := OldFit;
    Ind.FParameters := OldP;
   end;
  Finalize(OldP);
  Inc(FStatistics.Calculations);
  FStatistics.NewBest := FStatistics.NewBest or (i=1);
 end;
end;


procedure TGeneticAlgorithm.MakeOffsprings;
var
 i,P1i,P2i,ind,Ci: integer;
 P1,P2,C: TIndividual;
begin
 SetSelectionParametrs;
 try
  for i := 1 to FChildrensCount do
  begin
   if i = 1 then P1i := 0 else       //������ ������ ���� ���������
   P1i := random(FParentsCount);
   P2i := random(FParentsCount-1);

   if P2i >= P1i then Inc(P2i);
   if P1i > P2i then begin ind := P1i; P1i := P2i; P2i := ind; end;
   P1 := TIndividual(FPopulation[P1i]);
   P2 := TIndividual(FPopulation[P2i]);
   Ci := FPopulationSize-1 - FChildrensCount + i;
   C := TIndividual(FPopulation[Ci]);

   P1.MakeSexChromosome(C.FChrLg,C.FChrLd,CrossCountMean,MutationRate,TranslocationProbability,InversionProbability);
   P2.MakeSexChromosome(C.FChrRg,C.FChrRd,CrossCountMean,MutationRate,TranslocationProbability,InversionProbability);
   C.DecodePhenotype;

   C.FAdge := 0;

   with FStatistics  do
   begin
    Inc(Calculations);
    CrossoversAv := CrossoversAv + (P1.CrossoversCount+P2.CrossoversCount)*FStatCoeff;
    TranslocationsAv := TranslocationsAv + (P1.TranslocationsCount+P1.TranslocationsCount)*FStatCoeff;
    InversionsAv := InversionsAv + (P1.InversionsCount+P1.InversionsCount)*FStatCoeff;
    MutationsAv := MutationsAv + (P1.MutationsCount+P2.MutationsCount)*FStatCoeff;
   end;

   if UseEvaluationThreads
   then EvaluateIndividualInThread(C)
   else EvaluateIndividual(C);
  end;
 finally
  if UseEvaluationThreads
  then WaitForEvaluationThreads;
 end;
end;



procedure TGeneticAlgorithm.Sort;
  function Better(Ind: TIndividual; MidF: TData; MidC: boolean; MidA: integer; TType: TMaxMin): boolean; inline;
  begin
   with Ind do
   begin
    result := Condition and (Adge >= 0);
    if result and MidC and (MidA >= 0) then
    begin
     if TType = tcMinimize
     then result := Fitness < MidF
     else result := Fitness > MidF;
    end;
   end;
  end;

  function Worse(Ind: TIndividual; MidF: TData; MidC: boolean; MidA: integer; TType: TMaxMin): boolean; inline;
  begin
   with Ind do
   begin
    result := MidC and (MidA >= 0);
    if result and Condition and (Adge >= 0) then
    begin
     if TType = tcMinimize
     then result := Fitness > MidF
     else result := Fitness < MidF;
    end;
   end;
  end;

  procedure QuickSort(iLo, iHi: integer);
  var Lo, Hi: integer;
      MidF: TData;
      MidC: boolean;
      MidA: integer;
  begin
    Lo := iLo;
    Hi := iHi;
    with TIndividual(Population[(Lo + Hi) div 2]) do
    begin
     MidF := Fitness;
     MidC := Condition;
     MidA := Adge;
    end;
    repeat
      while Better(TIndividual(FPopulation[Lo]),MidF,MidC,MidA,MaxMinTask) do Inc(Lo);
      while Worse(TIndividual(FPopulation[Hi]),MidF,MidC,MidA,MaxMinTask) do Dec(Hi);
      if Lo <= Hi then
      begin
       if Lo <> Hi then
       FPopulation.Exchange(Lo,Hi);
       Inc(Lo);
       Dec(Hi);
      end;
    until Lo > Hi;
    if Hi > iLo then QuickSort(iLo, Hi);
    if Lo < iHi then QuickSort(Lo, iHi);
  end;

begin
 QuickSort(0, PopulationSize-1);

 with FStatistics do
 NewBest := NewBest or (BestIndividual <> Population[0]);
 BestIndividual := TIndividual(Population[0]);
end;


function DaethProbability(a,MLL: integer): single; inline;
begin
 if a = 0
 then result := 0
 else result := 0.01*exp(8*ln(a/MLL)); //������� ��������
 if result > 1 then result := 1;
end;


procedure TGeneticAlgorithm.AgeControl;
var i: integer;
begin
 FDeadCount := 0;
 for i := 0 to FPopulationSize-1 do
 with TIndividual(FPopulation[i]) do
 begin
  Inc(FAdge);
  if (DeathCondition = dcDeath) and
     (random < DaethProbability(FAdge,MeanLifeLengh)) and
     (FPopulationSize > 2) then
  begin
   FAdge := -1;
   Inc(FDeadCount);
  end;
 end;
end;



procedure TGeneticAlgorithm.DoNextGeneration;
var
 StartTime: TDateTime;
const
 WinFreezTime = 5/(24*60*60); // ~����� ��������� ��������� ��� ������ ����������
begin
 StartTime := Now;
 if PopulationValid
 then
  begin
   FStatistics.NewBest := false;
   if ReevaluateCondition = ecNow
   then
    begin
     EvaluateAll;
     Sort;
     ReevaluateCondition := ecNew;
    end
   else
    begin
     ClearGenerationStatistics;
     AgeControl;
     Sort;
     MakeOffsprings;
     if ReevaluateCondition = ecAll then EvaluateAll;
     if OntogenesisCondition = ocOntAll then OntAll;
     Inc(FStatistics.Epochs);
     CalcStatistics;
    end;

   with FStatistics do
   begin
    EpochTime := Now - StartTime;
    if HibernationStartTime > 0 then            //���� ������� � ����������
    begin
     if HibernationStopTime > 0
     then                                        //��������� ���������
      begin
       EpochTime := EpochTime - (HibernationStopTime-HibernationStartTime);
       HibernationStartTime := 0;
       HibernationStopTime := 0;
      end
     else                                        //���-�� �� ��
      begin
       if Now - HibernationStartTime < WinFreezTime
       then                                       //��� �� �����
        HibernationStartTime := Now
       else                                       //��������� ������
        begin
         if Epochs > 1
         then EpochTime := EvolutionTime/(Epochs-1)
         else EpochTime := HibernationStartTime - StartTime;
         HibernationStartTime := 0;
        end
      end;
    end;
    if EpochTime > 0 then                      //������ �� �������� ����� �����
    EvolutionTime := EvolutionTime + EpochTime;
   end;
  end
 else
  raiseError('Population is not full');
end;


procedure TGeneticAlgorithm.CalcChromosomeParams;
var
 i: Integer;
begin
 FChromosomeLen := 0;
 for i := 0 to FParametersVectorLen - 1 do
 begin
  FLinkGroupShift[i] := FChromosomeLen;
  FChromosomeLen := FChromosomeLen + FChromosomeMap[i];
 end;
end;


procedure TGeneticAlgorithm.ClearGenerationStatistics;
begin
 with FStatistics do
 begin
  CrossoversAv := 0;
  MutationsAv := 0;
  InversionsAv := 0;
  TranslocationsAv := 0;
 end;
end;


procedure TGeneticAlgorithm.CalcStatistics;
var i: integer;
begin
 with FStatistics do
 begin
  ParentsCount :=  FParentsCount;
  ChildrensCount := FChildrensCount;
  DeadCount := FDeadCount;

  EfRadius := 0;
  MeanAge := TIndividual(Population[0]).Adge;
  if ParentsCount > 1 then
  begin
   i := 1;
   repeat
    EfRadius := EfRadius + Distance(TIndividual(Population[0]).parameters,TIndividual(Population[i]).parameters);
    MeanAge := MeanAge + TIndividual(Population[i]).Adge;
    inc(i);
   until i >= ParentsCount;
   EfRadius := EfRadius/(ParentsCount-1);
   MeanAge := MeanAge/ParentsCount;
  end;

  DifFitn := abs(TIndividual(Population[0]).Fitness-TIndividual(Population[ParentsCount-1]).Fitness);
 end;
end;


initialization
 Randomize;
 FCSection := TCriticalSection.Create;
 GAProcessorsCount := 2;
 GAThreadsPerProcessor := 3;

finalization
 FCSection.Free;

END.


