{***********************************}
{*                                 *}
{* ћодуль базовых типов и процедур *}
{*  дл€ математических вычислений  *}
{*                                 *}
{***********************************}
{
20.03.2008
 + функции копировани€ матриц и векторов
20.03.2008
23.02.2006
17.10.2002
}

unit BaseMath;

interface
uses classes;
Type
 TData = double;
 TVector = array of TData;
 TMatrix = array of TVector;
 TIntVector = array of integer;

 function CopyVector(const A: TVector):TVector; inline;
 function CopyMatrix(const A: TMatrix):TMatrix; inline;
 procedure SaveTVectorToStream(Stream: TStream; var Vector: TVector); inline;
 procedure LoadTVectorFromStream(Stream: TStream; var Vector: TVector); inline;
 procedure SaveTIntVectorToStream(Stream: TStream; var IntVector: TIntVector); inline;
 procedure LoadTIntVectorFromStream(Stream: TStream; var IntVector: TIntVector); inline;

 procedure DeleteRowFromMatrixRO(var M: TMatrix; pos: integer); inline;
 procedure DeleteFromVectorRO(var V: TVector; pos: integer); overload;  inline;
 procedure DeleteFromVectorRO(var V: TIntVector; pos: integer); overload; inline;

implementation

function CopyVector(const A: TVector):TVector;
begin
 Result := copy(A);
end;


function CopyMatrix(const A: TMatrix):TMatrix;
var
 i: integer;
begin
 SetLength(Result, High(A)+1);
 for i := 0 to high(A) do
 begin
  Result[i] := copy(A[i]);
 end;
end;


procedure SaveTVectorToStream(Stream: TStream; var Vector: TVector);
var l: integer;
begin
 l := length(Vector);
 Stream.WriteBuffer(l, SizeOf(integer));
 Stream.WriteBuffer(Vector[0], SizeOf(TData)*l);
end;


procedure LoadTVectorFromStream(Stream: TStream; var Vector: TVector);
var l: integer;
begin
 Stream.ReadBuffer(l, SizeOf(integer));
 SetLength(Vector,l);
 Stream.ReadBuffer(Vector[0], SizeOf(TData)*l);
end;


procedure SaveTIntVectorToStream(Stream: TStream; var IntVector: TIntVector);
var l: integer;
begin
 l := length(IntVector);
 Stream.WriteBuffer(l, SizeOf(integer));
 Stream.WriteBuffer(IntVector[0], SizeOf(integer)*l);
end;


procedure LoadTIntVectorFromStream(Stream: TStream; var IntVector: TIntVector);
var l: integer;
begin
 Stream.ReadBuffer(l, SizeOf(integer));
 SetLength(IntVector,l);
 Stream.ReadBuffer(IntVector[0], SizeOf(integer)*l);
end;


procedure DeleteRowFromMatrixRO(var M: TMatrix; pos: integer);
//!!пор€док следовани€ нарушаетс€!!
begin
 Finalize(M[pos]);
 M[pos] := M[High(M)];
 SetLength(M,Length(M)-1);
end;


procedure DeleteFromVectorRO(var V: TVector; pos: integer); overload;
//!!пор€док следовани€ нарушаетс€!!
begin
 V[pos] := V[High(V)];
 SetLength(V,Length(V)-1);
end;


procedure DeleteFromVectorRO(var V: TIntVector; pos: integer); overload;
//!!пор€док следовани€ нарушаетс€!!
begin
 V[pos] := V[High(V)];
 SetLength(V,Length(V)-1);
end;


end.



