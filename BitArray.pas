{*********************************************}
{*                                           *}
{*  Модуль динамического битового  массива   *}
{*                                           *}
{*********************************************}
{* Based on Delphi Visual Component Library  *}
{* Classes Unit                              *}
{*********************************************}
{* первый индекс в векторе - 1               *}
{*********************************************}
{
27.03.16
 * восстановлена совместимость с x64
31.10.09
 - серьезная ошибка в CopyBitsTo, приводившая при определенных
   сочетаниях параметров к тому, что данные не копировались
03.10.09
 * небольшая оптимизация кода в CopyBitsTo
21.03.09
 * небольшая оптимизация кода в CopyBitsTo
06.01.08
 * введение метода CopyBitsTo для ускорения операций с TBitVector
18.10.02;
30.01.99;
30.01.98;
}

unit BitArray;
{$DEFINE X86ASM}

interface

uses Classes,SysUtils,Math;

Type
  ERangeInvalid = class(Exception);

  TBitVector = class
  private
   FSize: integer;
   FMemSize: integer;
   FBits: Pointer;
   procedure Error;  {$IF CompilerVersion >= 18.0} inline; {$IFEND}
   procedure SetSize(Value: integer);
   procedure SetBit(Index: integer; Value: Boolean);
   function GetBit(Index: integer): Boolean;
  public
   property Bits[Index: integer]: Boolean read GetBit write SetBit; default;
   property Size: integer read FSize write SetSize;
   procedure CopyBitVectorFrom(Source: TBitVector);  {$IF CompilerVersion >= 18.0} inline; {$IFEND}
   procedure CopyBitsTo(Dest: TBitVector; IndexSrc, Count, IndexDst: integer);
   procedure SaveToStream(Stream: TStream);
   constructor LoadFromStream(Stream: TStream);
   destructor Destroy; override;
  end;

implementation

{-----------------------------------TBitVector---------------------------------}
const
  BitsPerInt = SizeOf(Integer) * 8;


destructor TBitVector.Destroy;
begin
  SetSize(0);
  inherited Destroy;
end;


procedure TBitVector.Error;
begin
  raise ERangeInvalid.Create('Bits index out of range');
end;


function Min(X, Y: integer): integer; {$IF CompilerVersion >= 18.0} inline; {$IFEND}
begin
  Result := X;
  if X > Y then Result := Y;
end;


procedure TBitVector.SetSize(Value: integer);
var
  NewMem: Pointer;
  NewMemSize: integer;
begin
  if Value <> Size then
  begin
    if Value < 0 then Error;
    NewMemSize := ((Value + BitsPerInt - 1) div BitsPerInt) * SizeOf(integer);
    if NewMemSize <> FMemSize then
    begin
      NewMem := nil;
      if NewMemSize <> 0 then
      begin
        GetMem(NewMem, NewMemSize);
        FillChar(NewMem^, NewMemSize, 0);
      end;
      if FMemSize <> 0 then
      begin
        if NewMem <> nil then
          Move(FBits^, NewMem^, Min(FMemSize, NewMemSize));
        FreeMem(FBits, FMemSize);
      end;
      FBits := NewMem;
      FMemSize := NewMemSize;
    end;
    FSize := Value;
  end;
end;


procedure TBitVector.SetBit(Index: integer; Value: Boolean); assembler;
{$IFNDEF X86ASM}
var
  LRelInt: PInteger;
  LMask: Integer;
begin
  dec(Index);
  if (Index >= FSize) or (Index < 0) then
    Error;
  { Calculate the address of the related integer }
  { Calculate the address of the related integer }
  LRelInt := FBits;
  Inc(LRelInt, Index div BitsPerInt);

  { Generate the mask }
  LMask := (1 shl (Index mod BitsPerInt));

  { Update the integer }
  if Value then
    LRelInt^ := LRelInt^ or LMask
  else
    LRelInt^ := LRelInt^ and not LMask;
end;
{$ELSE X86ASM}
asm
        DEC     Index
        CMP     Index,[EAX].FSize
        JAE     TBitVector.Error

@@1:    MOV     EAX,[EAX].FBits
        OR      Value,Value
        JZ      @@2
        BTS     [EAX],Index
        RET

@@2:    BTR     [EAX],Index
        RET
end;
{$ENDIF X86ASM}


function TBitVector.GetBit(Index: integer): Boolean; assembler;
{$IFNDEF X86ASM}
var
  LRelInt: PInteger;
  LMask: Integer;
begin
  dec(Index);
  if (Index >= FSize) or (Index < 0) then
    Error;
  { Calculate the address of the related integer }
  LRelInt := FBits;
  Inc(LRelInt, Index div BitsPerInt);

  { Generate the mask }
  LMask := (1 shl (Index mod BitsPerInt));
  Result := (LRelInt^ and LMask) <> 0;
end;
{$ELSE X86ASM}
asm
        DEC     Index
        CMP     Index,[EAX].FSize
        JAE     TBitVector.Error
        MOV     EAX,[EAX].FBits
        BT      [EAX],Index
        SBB     EAX,EAX
        AND     EAX,1
end;
{$ENDIF X86ASM}


procedure TBitVector.CopyBitVectorFrom(Source: TBitVector);
begin
 Move(Source.FBits^,FBits^,FMemSize);
end;


procedure TBitVector.CopyBitsTo(Dest: TBitVector; IndexSrc, Count, IndexDst: integer);
const
 bsz = sizeof(integer);
var
 n,sbt,dbt,n_Count: integer;
 sdw,ddw,hdw,
 i,buf,mask1,mask2,m: cardinal;
 psdw,pddw: pointer;
begin
 sdw := (IndexSrc-1) div BitsPerInt;
 sbt := (IndexSrc-1) mod BitsPerInt;
 ddw := (IndexDst-1) div BitsPerInt;
 dbt := (IndexDst-1) mod BitsPerInt;
 hdw := (Count - 1) div BitsPerInt;

 m := 1 shl dbt;      //00010000
 mask1 := m-1;        //00001111
 mask2 := not mask1;  //11110000

 n := 0;
 i := 0;
 repeat //for i := 0 to hdw do
  inc(n,BitsPerInt);
  n_Count := n-Count;

  //чтение источника
  psdw := pointer(cardinal(FBits)+(sdw+i)*bsz);
  buf := cardinal(psdw^) shr sbt;
  if (sbt > 0) and (n_Count < sbt) then //...n-sbt < Count
  begin
   psdw := pointer(cardinal(FBits)+(sdw+i+1)*bsz);
   buf := buf or (cardinal(psdw^) shl (BitsPerInt-sbt));
  end;

  //маски для записи получателя
  if n_Count > 0 then //n > Count ->  конец, прочитано больше, чем нужно
  begin
   m := 1 shl (BitsPerInt - n_Count);
   buf := buf and (m-1); //00001111
   if n_Count > dbt
   then
    begin
     m := 1 shl (dbt + BitsPerInt - n_Count);
     mask1 := mask1 or not(m-1);  //11110000
    end
   else
    if n_Count < dbt
    then
     begin
      m := 1 shl (dbt - n_Count);
      mask2 := mask2 or not(m-1);  //11110000
     end;
  end;

  //запись получателя
  pddw := pointer(cardinal(Dest.FBits)+(ddw+i)*bsz);
  cardinal(pddw^) := (cardinal(pddw^) and mask1) or (buf shl dbt);
  if (dbt > 0) and (n_Count < dbt)  then //...n-dbt < Count
  begin
   pddw := pointer(cardinal(Dest.FBits)+(ddw+i+1)*bsz);
   cardinal(pddw^) := (cardinal(pddw^) and mask2) or (buf shr (BitsPerInt-dbt));
  end;
  inc(i)
 until i > hdw;
end;


procedure TBitVector.SaveToStream(Stream: TStream);
begin
 Stream.WriteBuffer(FSize,sizeof(FSize));
 Stream.WriteBuffer(FBits^,FMemSize);
end;


constructor TBitVector.LoadFromStream(Stream: TStream);
var s: integer;
begin
 inherited Create;
 try
  Stream.ReadBuffer(s,sizeof(s));
  SetSize(s);
  Stream.ReadBuffer(FBits^,FMemSize);
 except
 on E: EOutOfMemory do
  begin
   E.Message := 'Can not load dynamic array';
   raise;
  end;
 end;
end;

END.
