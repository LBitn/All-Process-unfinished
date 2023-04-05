# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
module F2023
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
using HDF5
using Dates
using DelimitedFiles
using JLD
using Plots
using Statistics
using StatsBase
using DSP
using BinningAnalysis
using HistogramThresholding
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
export HDF5Content
export FindKeyRegex
export GetVarsHDF5
export ChunkSizeSpace
export OneSegment
export Digital2Analogue
export UniqueCount
export RemoveInfs
export Zplot
export ms2Frs
export stdΔV
export Neighbors
export convgauss
export RemoveSaturationExtrema
export FindContent
export MatrixFilter
export FindAllThrs
export SigmaData
export Energy
export Donoho
export OutlayerME
export LoadDict
export AllMs
export GetIndexes
export PFT_mt
export CSD
export ComponentesSP
export vecindad8
export CentrosMasa
export CMs
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
GaussianKernel = [
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
    0.00002292	0.00078634	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
    0.00038771	0.01330373	0.11098164	0.22508352	0.11098164	0.01330373	0.00038771
    0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
    0.00002292	0.00078633	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
];
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
# El operador de Laplace-Lindenberg
Laplacian01 = [ [ 0 1 0 ]; [ 1 -4 1 ]; [ 0 1 0 ] ];
Laplacian02 = [ [ 0.5 0 0.5 ]; [ 0 -2 0 ]; [ 0.5 0 0.5 ] ];
LaplacianKernel = ( ( 1 - ( 1 / 3 ) ) * Laplacian01 ) + ( ( 1 / 3 ) * Laplacian02 );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #

# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    GetGroupsHDF5( BRW::HDF5.File, g::String ) -> GroupsN::Vector{String}
        Extract the groups form a BRW open file.
"""
function GetGroupsHDF5( BRW::HDF5.File, g::String )
    GroupsN = [ ];
    try
        GroupsN = string.( g, "/", keys( BRW[ g ] ) );
    catch e
    end
    return GroupsN
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    GetAttrHDF5( BRW::HDF5.File, g::String ) -> AttrN::Vector{String}
        Extract the attributes form a BRW open file.
        using HDF5
"""
function GetAttrHDF5( BRW::HDF5.File, g::String )
    AttrN = [ ];
    aux = attributes( BRW[ g ] );
    try
        AttrN = string.( g, "/", keys( aux ) );
    catch e
    end
    return AttrN
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    HDF5Content( BRW::HDF5.File ) -> AllGroups::Vector{Any}
        Generates a list with the entire content of an HDF5 file.
        using GetGroupsHDF5, GetAttrHDF5
"""
function HDF5Content( BRW::HDF5.File )
    Groups = keys( BRW );
    AllGroups = [ ];
    AllAtributes = keys( attributes( BRW ) );
    while !isempty( Groups )
        GROUPS = [ ];
        ATTR = [ ];
        for g in Groups
            if typeof( BRW[ g ] ) == HDF5.Dataset
                push!( AllGroups, g )
            else
                push!( GROUPS, GetGroupsHDF5( BRW, g ) );
                AllAtributes = vcat( AllAtributes, GetAttrHDF5( BRW, g ) );
            end
        end
        Groups = vcat( GROUPS... );
        push!( AllGroups, Groups );
    end
    AllGroups = unique( vcat( AllGroups... ) );
    AllAtributes = unique( AllAtributes );
    return AllGroups, AllAtributes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ExtractRawDataset( BRW::HDF5.File, AllGroups::Vector{Any} ) -> Raw::String, NoRaw::String
        Detects which is the largest dataset and asigns it to Raw data.
"""
function ExtractRawDataset( BRW::HDF5.File, AllGroups::Vector{Any} )
    Types = Vector{ String }( undef, length( AllGroups ) );
    for g = 1:length( AllGroups )
        Types[ g ] = string( typeof( BRW[ AllGroups[ g ] ] ) );
    end
    AllDataSets = AllGroups[ Types .== "HDF5.Dataset" ];
    aux = zeros( Int, length( AllDataSets ) );
    for i = 1:length( AllDataSets )
        aux[ i ] = length( BRW[ AllDataSets[ i ] ] );
    end
    Raw = AllDataSets[ aux .== maximum( aux ) ][ 1 ];
    NoRaw = AllDataSets[ aux .!= maximum( aux ) ];
    aux = aux[ aux .!= maximum( aux ) ];
    NoRaw = NoRaw[ aux .!= 0 ];
    return Raw, NoRaw
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ExtractValues( BRW::HDF5.File,  AllAtributes::Vector{String}, NoRaw::Vector{Any} ) -> D::Dict{Any, Any}
        Extracts the values from the BRW file, on Float or Int format if posible.
        using HDF5
"""
function ExtractValues( BRW::HDF5.File,  AllAtributes::Vector{String}, NoRaw::Vector{Any} )
    D = Dict( ); e = [ ];
    for g in NoRaw
        try
            D[ g ] = Float64.( read( BRW[ g ] ) );
        catch e
            D[ g ] = read( BRW[ g ] );
        end
    end
    for g in keys( D )
        if length( D[ g ] ) .== 1
            D[ g ] = D[ g ][ 1 ];
        end
        try
            D[ g ] = Int64( D[ g ] );
        catch e
        end
    end
    for g in AllAtributes
        aux00 = basename( g ); aux01 = dirname( g );
        if isempty( aux01 )
            D[ g ] = read_attribute( BRW, aux00 );
        else
            D[ g ] = read_attribute( BRW[ aux01 ], aux00 );
        end
    end
    return D
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    FindKeyRegex( key::String, D::Dict{Any, Any} ) -> Vector{String}
        Searchs the entries of the Dictionary D for the key word using Word Match.
"""
function FindKeyRegex( key::String, D::Dict{Any, Any} )
    key = Regex( ".+$key.+" );
    ok = keys( D );
    aux = match.( key, ok );
    aux01 = aux[ aux .!= nothing ];
    okok = [ ];
    if !isempty( aux01 )
        for i in aux01
            push!( okok, i.match );
        end
        return okok
    else
        println( "There is no match for that keyword into the Dicctionary" );
    end
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    FindParameters( D::Dict{Any, Any} ) -> D::Dict{Any, Any}
        Search for the keywords "Ch" and "Std" to extract the values or generate them if its missing from the original file.
        using FindKeyRegex
"""
function FindParameters( D::Dict{Any, Any} )
    aux0 = FindKeyRegex( "Ch", D );
    if !isempty( aux0 )
        aux1 = [ ];
        for i in aux0
            try
                push!( aux1, size( D[ i ], 1 ) );
            catch e
                push!( aux1, 0 );
            end
        end
        nChs = size( D[ aux0[ aux1 .== maximum( aux1 ) ][ ] ], 1 );
    else
        nChs = 4096;
    end
    aux0 = FindKeyRegex( "Std", D );
    if !isempty( aux0 )
        aux1 = [ ];
        for i in aux0
            try
                push!( aux1, size( D[ i ], 1 ) );
            catch e
                push!( aux1, 0 );
            end
        end
        STD = D[ aux0[ aux1 .== maximum( aux1 ) ][ ] ];
    else
        STD = [ ];
    end
    D[ "nChs" ] = nChs;
    D[ "STD" ] = STD;
    return D
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ExpDate2Str( Variables::Dict ) -> Variables::Dict
        Extracting the date of the BRW creation
        using Dates
"""
function ExpDate2Str( Variables::Dict )
    X = Variables[ "ExperimentDateTime" ];
    Dt = split( X, ":" );
    Dt[ end ] = string( round( Int, parse( Float64, replace(
        Dt[ end ], r"[A-Z]" => "" ) ) ) );
    newDt = String( "" );
    for i in Dt
        newDt = string( newDt, ":", i );
    end
    newDt = newDt[ 2:end ];
    X = Dates.DateTime( newDt );
    Variables[ "ExperimentDateTime" ] = string( X );
    T = Dates.format( X, RFC1123Format );
    println( "Creation Date: ", T )
    return Variables
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ExpSett2Dict( Variables::Dict ) -> Variables::Dict
        Extracting the contents of the ExperimentSettings dictionary
"""
function ExpSett2Dict( Variables::Dict )
    ExperimentSettings = Variables[ "ExperimentSettings" ][ 1 ];
    t = split( ExperimentSettings,"\r\n" );
    t = replace.( t, "  " => "", "{" => "", "}" => "", '"' => "" );
    x = [ ];
    for i = 1:length( t )
        if !isempty( collect( eachmatch( r"[a-z]", t[ i ] ) ) )
            push!( x, true );
        else
            push!( x, false );
        end
    end
    t = t[ Bool.( x ) ]; t = split.( t, ": " );
    D = Dict( );
    for i in t
        if !( i[ 2 ] == "" )
            aux = i[ 2 ];
            try
                aux = replace( i[ 2 ], "," => " ", "[" => "", "[]" => "","]" => "", " " => "" )
            catch e
            end
            if ( aux != "" ) && ( aux != " " )
                aux = aux;
                try
                    aux = parse( Float64, aux );
                catch
                    aux = replace( aux, " " => "" );
                end
                D[ i[ 1 ] ] = aux;
            end
        end
    end
    delete!( Variables, "ExperimentSettings" );
    Variables = merge( Variables, D );
    return Variables
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    CleanDictionary( D::Dict{Any, Any} ) -> D::Dict{Any, Any}
        Clean the entry dictionary form empty or null values.
"""
function CleanDictionary( D::Dict{Any, Any} )
    X = String.( keys( D ) )[ values( D ) .== "null" ];
    for x in X
        delete!( D, x );
    end
    X = [ ];
    for k in keys( D )
        try
            if isempty( D[ k ] )
                push!( X, k );
            end
        catch e
        end
    end
    for x in X
        delete!( D, x );
    end
    return D
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Vars2TXT( D::Dict{Any, Any}, lim = 100 ) -> TXT::Array{String}
        Transforms the D dictionary entries to a text array with limit word count to 100 or lim
"""
function Vars2TXT( D::Dict{Any, Any}, lim = 100 )
    x = [ ];
    for i in values( D )
        if length( i ) <= lim
            push!( x, 1 )
        else
            push!( x, 0 )
        end
    end
    NK = string.( keys( D ) )[ Bool.( x ) ];
    TXT = Dict( ); for i in NK; TXT[ i ] = D[ i ]; end
    return TXT
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    PATHSFILES( FILEBRW::String, Variables::Dict ) -> FILEPATHS::String, FILEVARS::String, PATHS::Dict{String, String}
        Generate the folders, and save the variable file and pats file. Also displays de brw description.
        using DelimitedFiles, JLD
        using Vars2TXT
"""
function PATHSFILES( FILEBRW::String, Variables::Dict )
    PATHBRWs = dirname( FILEBRW );
    PATHMAIN = joinpath( dirname( PATHBRWs ), split( basename( FILEBRW ), "." )[ 1 ] );
    PATHINFO = joinpath( PATHMAIN, "Info" ); mkpath( PATHINFO );
    FILEVARS = joinpath( PATHINFO, "Variables.jld" );
    FILEVARSTXT = joinpath( PATHINFO, "Variables.txt" );
    FILEPATHS = joinpath( PATHINFO, "Paths.jld" );
    PATHS = Dict(
        "PATHMAIN" => PATHMAIN,
        "PATHINFO" => PATHINFO,
        "PATHBRWs" => PATHBRWs
        );
    TXT = Vars2TXT( Variables );
    writedlm( FILEVARSTXT, TXT );
    save( FILEVARS, "Variables", Variables );
    save( FILEPATHS, "PATHS", PATHS );
    BRWsize = Variables[ "BRWsizeGB" ];
    BRWname = basename( Variables[ "BRWname" ] );
    println( "You are now working on the new main path: ", PATHMAIN );
    println( "With the file: ")
    print( BRWname, " : ", Variables[ "Description" ] );
    println( " HDF5 file size: $BRWsize GB" );
    return FILEPATHS, FILEVARS, PATHS
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    GetVarsHDF5( FILEBRW::String ) -> Variables::Dict, FILEPATHS::String, FILEVARS::String
        Extracts all contents from a HDF5 file of all versions, except the Raw dataset, which only shows the string with the location
        using HDF5, Dates, DelimitedFiles, JLD
        using HDF5Content, ExtractRawDataset, ExtractValues, FindParameters, ExpSett2Dict, ExpDate2Str, CleanDictionary, PATHSFILES
"""
function GetVarsHDF5( FILEBRW::String )
    BRW = h5open( FILEBRW, "r" );
    AllGroups, AllAtributes = HDF5Content( BRW );
    Raw, NoRaw = ExtractRawDataset( BRW, AllGroups );
    Variables = ExtractValues( BRW,  AllAtributes, NoRaw );
    Variables = FindParameters( Variables );
    dset = BRW[ Raw ];
    BRWsize =  ( stat( FILEBRW ).size ) / ( 1024 * 1024 * 1024 )
    DataSetSize = (
        sizeof( typeof( dset[ 1 ] ) ) * size( dset )[ 1 ] ) / ( 1024 * 1024 * 1024 );
    Variables[ "Raw" ] = Raw;
    Variables[ "BRWname" ] = BRW.filename;
    Variables[ "BRWsizeGB" ] = BRWsize;
    Variables[ "DataSetSize" ] = DataSetSize;
    NewVars = Dict( ); K = keys( Variables );
    for k in K; NewVars[ basename( k ) ] = Variables[ k ]; end; Variables = NewVars;
    if ( "ExperimentSettings" in keys( Variables) )
        Variables = ExpSett2Dict( Variables);
        Variables = ExpDate2Str( Variables);
    end
    Variables = CleanDictionary( Variables );
    FILEPATHS, FILEVARS, PATHS = PATHSFILES( FILEBRW, Variables );
    cd( PATHS[ "PATHMAIN" ] )
    close( BRW )
    return Variables, FILEPATHS, FILEVARS
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ChunkSizeSpace( Variables::Dict, limupper::Real ) -> σ::Int64
        Estabish the number of segmets to cut form each brw file considering the upper limit in GB size limupper for arrays in Float16 format
"""
function ChunkSizeSpace( Variables::Dict, limupper::Real )
    NRecFrames = Variables[ "NRecFrames" ];
    SamplingRate = Variables[ "SamplingRate" ];
    σ = 5
    finalsize = Variables[ "DataSetSize" ] / σ;
    while ( finalsize > limupper )
        σ = σ + 1
        finalsize = Variables[ "DataSetSize" ] / σ;
        while NRecFrames % σ != 0
            σ = σ + 1
        end
    end
    finalsize = Variables[ "DataSetSize" ] / σ;
    finaltime = ( NRecFrames / σ ) / SamplingRate; # sec
    fs = round( finalsize, digits = 3 );
    ft = round( finaltime, digits = 3 );
    println( "$σ segments of $fs GB and $ft seconds each" );
    return σ
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    OneSegment( Variables::Dict, n::Int64, nSegments::Int ) -> BIN::::Matrix{Float64}
        Takes just one segment of the BRW file.
        using HDF5
"""
function OneSegment( Variables::Dict, n::Int64, nSegments::Int )
    nChs = Variables[ "nChs" ];
    NRecFrames = Variables[ "NRecFrames" ];
    binlenght = Int( NRecFrames / nSegments );
    Raw = h5open( Variables[ "BRWname" ], "r" )[ Variables[ "Raw" ] ];
    BIN = Array{ UInt16 }( undef, nChs, binlenght );
    for frame = ( ( ( n - 1 ) * binlenght ) + 1 ): binlenght * n
        BIN[ :,
            ( frame - ( binlenght*( n - 1 ) ) ) ] .=
                Raw[ ( ( ( frame - 1 ) * nChs ) + 1 ): nChs * frame ];
    end
    return BIN
    close( Raw )
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::Matrix{UInt16} ) -> BIN::Matrix{Float64}
        Conversion of raw data extracted from the brw file to voltage values (μV) for Matrix format acording to the equation
        Voltage = ( RawData + ADCCountsToMV ) * MVOffset
"""
function Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::Matrix{UInt16} )
    SignalInversion = Variables[ "SignalInversion" ];
    MinVolt = Variables[ "MinVolt" ];
    MaxVolt = Variables[ "MaxVolt" ];
    BitDepth = Variables[ "BitDepth" ];
    MVOffset = SignalInversion * MinVolt;
    ADCCountsToMV = ( SignalInversion * ( MaxVolt - MinVolt ) ) / ( 2^BitDepth );
    BIN = @. MVOffset + ( DigitalValue * ADCCountsToMV );
    return BIN
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::UInt16 ) -> AnalogValue::Float64
    Conversion of raw data extracted from the brw file to voltage values (μV) for one entrie acording to the equation
        Voltage = ( RawData + ADCCountsToMV ) * MVOffset
"""
function Digital2Analogue( Variables::Dict{ Any, Any }, DigitalValue::UInt16 )
    SignalInversion = Variables[ "SignalInversion" ];
    MinVolt = Variables[ "MinVolt" ];
    MaxVolt = Variables[ "MaxVolt" ];
    BitDepth = Variables[ "BitDepth" ];
    MVOffset = SignalInversion * MinVolt;
    ADCCountsToMV = ( SignalInversion * ( MaxVolt - MinVolt ) ) / ( 2^BitDepth );
    AnalogValue = MVOffset + ( DigitalValue * ADCCountsToMV );
    return AnalogValue
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    UniqueCount( data::Array ) -> Count::Vec{Int64}
        Counts the number of unique values for each row of an array
"""
function UniqueCount( data::Array )
    X, Y = size( data );
    Count = Array{ Int64 }( undef, X );
    [ Count[ x ] = length( unique( round.( data[ x, : ], digits = 2 ) ) ) for x in 1:X ];
    return Count
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    RemoveInfs( data::VecOrMat )
        Remove the -Infs and Infs entries. Replace with maximum and minimum nonInf value
"""
function RemoveInfs( data::VecOrMat )
    m = minimum( data[ data .!= -Inf ] );
    M = maximum( data[ data .!= Inf ] );
    data[ data .== -Inf ] .= m;
    data[ data .== Inf ] .= M;
    return data
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Z0( X::VecOrMat, nChs::Int64 ) -> Z::Matrix{Int64}
        using Plots
"""
function Z0( X::VecOrMat, nChs::Int64 )
    X = Int.( vec( X ) );
    Z = zeros( Int, nChs );
    n = Int( sqrt( nChs ) );
    Z[ X ] .= Z[ X ] .+ 1;
    Z = reverse( reshape( Z, n, n )', dims = 1 );
    return Int.( Z )
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ZW( X::VecOrMat ) -> Z::Matrix{typeof(X)}
        using Plots
"""
function ZW( X::VecOrMat )
    X = vec( X );
    n = Int( sqrt( length( X ) ) );
    Z = reverse( reshape( X, n, n )', dims = 1 );
    return Z
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Zplot( Z::VecOrMat, which::String, cm = :greys ) -> F::Plot
        which = "0", entry vector of numbers for a Back and White plot
        using Z0, ZW
        using Plots
"""
function Zplot( Z::VecOrMat, which::String, cm = :greys, nChs = 4096 )
    if which == "0"
        Z = Z0( Z, nChs );
    elseif which == "W"
        Z = ZW( Z );
    end
    F = heatmap( Z, aspect_ratio = 1, c = cm, axis = ( [ ], false ), wsize = ( 400, 400 ) );
    if which == "0"
        plot!( F, cbar  = :none )
    end
    return F
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    ms2Frs( Variables::Dict, time::Real ) -> x::Float64
        For conversion from ms to an integer number of frames
"""
function ms2Frs( Variables::Dict, time::Real )
    SamplingRate = Variables[ "SamplingRate" ];
    if time != 0; x = ceil( Int, ( time * SamplingRate ) / 1000 ); else; x = 1; end
    return x
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    stdΔV( Variables::Dict, BIN::Matrix{Float64}, ΔT::Real ) -> STD::Vector{Float64}
        ΔT en msec
        using Statistics
        using ms2Frs
"""
function stdΔV( Variables::Dict, BIN::Matrix{Float64}, ΔT::Real )
    ΔT = ms2Frs( Variables, ΔT );
    STD = vec( std( ( BIN - circshift( BIN, ( 0, ΔT )  ) ), dims = 2 ) );
    return STD
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Neighbors( C::Int64, d::Int64 ) -> A::Array{ Int64 }, v::Vector{ Int64 }
        A = Array( ( d*2 ) + 1, ( d * 2 ) + 1 ),
        v = vec( 2*( ( d * 2 ) + 1 ) - 1 );
        The d-neighborhood is calculated from the channel ( C ) as a center
        A = array where C is the center and is in chip order
        v = same neighboring channels as A but in vector form and without C ( 8 channels for d = 1 )
"""
function Neighbors( center::Int64, d::Int64 )
    Layout = reverse( reshape( collect( 1:4096 ), 64, 64 )', dims = 1 );
    x = findall( Layout .== center )[ ][ 2 ];
    y = findall( Layout .== center )[ ][ 1 ];
    aux = [ ( x - d ),( x + d ), ( y - d ), ( y + d ) ]
    aux[ aux .< 1 ] .= 1; aux[ aux .> 64 ] .= 64;
    A = Layout[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    v = vec( A )[ vec( A ) .!= center ];
    return A, v
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    convgauss( sigma::Real, h::Vector ) -> hg::Vector'
        Convolution between a gaussian function and a vector
"""
function convgauss( sigma::Real, h::Vector )
    cte = ( 1 / ( sigma * sqrt( 2 * pi ) ) );
    potencia( x ) = -( ( x ^ 2 ) / ( 2 * ( sigma ^ 2 ) ) );
    gx( x ) = cte * exp( potencia( x ) );
    g = gx.( h );
    hg = h .* g;
    return hg
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    RemoveSaturationExtrema( Constants::Dict, BIN::Matrix ) -> BIN::Matrix{Float64}, Discarded::Dict{Any,Any}
        Detects voltage values above the maximum limit and below the minimum limit permited.
        Saturations above the permited percetage, all channel is reconstructed with the 8-neighborhood mean, if is below percetage, is filled with random values from the same channel
        using StatsBase
        using Neighbors
"""
function RemoveSaturationExtrema( Constants::Dict, BIN::Matrix )
    # ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~  #
    maxvolt = Constants[ "maxvolt" ];
    minvolt = Constants[ "minvolt" ];
    limite_saturacion = Constants[ "limite_saturacion" ];
    minchannels = Constants[ "minchannels" ];
    maxradius = Constants[ "maxradius" ];
    nChs = size( BIN, 2 );
    channels = collect( 1:nChs );
    Discarded = Dict( );
    # ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~  #
    aux = findall( BIN .>= maxvolt .|| BIN .<= minvolt );
    SatChannels = getindex.( aux, [ 1 ] );
    SatFrames = getindex.( aux, [ 2 ] );
    aux01 = unique( SatChannels );
    AllFrames = [ ];
    AllChannels = [ ];
    for ch in aux01
        push!( AllChannels, ch );
        push!( AllFrames, SatFrames[ SatChannels .== ch ] );
    end
    AllFrames = AllFrames[ sortperm( AllChannels ) ];
    AllChannels = sort( AllChannels );
    ValidChs = setdiff( channels, AllChannels );
    SatLevel = length.( AllFrames ) / size( BIN, 2 );
    Empties = AllChannels[ SatLevel .>= limite_saturacion ];
    if isempty( Discarded )
        Discarded[ "Channels" ] = AllChannels;
        Discarded[ "Frames" ] = AllFrames;
        Discarded[ "Empties" ] = Empties;
    end
    # ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~  #
    if !isempty( Empties )
        for emptie in Empties
            radius = 1;
            _, neighs = Neighbors( emptie, radius );
            ValidNeighbors = [ ];
            ValidNeighbors = intersect( neighs, ValidChs );
            A = ( length( ValidNeighbors ) .>= minchannels );
            while !A && radius <= maxradius
                radius = radius + 1;
                _, neighs = Neighbors( emptie, radius );
                ValidNeighbors = intersect( neighs, ValidChs );
                A = ( length( ValidNeighbors ) .>= minchannels );
            end
            BIN[ emptie, : ] = mean( BIN[ ValidNeighbors, : ], dims = 1 );
        end
    end
    # ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~  #
    aux = AllChannels .∉ [ Empties ];
    SatFrames = AllFrames[ aux ];
    SatChannels = AllChannels[ aux ];
    for i = 1:length( SatChannels );
        satch = SatChannels[ i ];
        satfrs = SatFrames[ i ];
        Nsatfrs = length( satfrs );
        ValidFrs = setdiff( 1:size( BIN, 2 ), satfrs );
        BIN[ satch, satfrs ] = sample( BIN[ satch, ValidFrs ], Nsatfrs );
    end
    # ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~  #
    return BIN, Discarded
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    FindContent( thing::String, FILEBRW::String = "", WORKDIR::String = "" )
"""
function FindContent( thing::String, FILEBRW::String = "", WORKDIR::String = "" )
    if FILEBRW != ""
        WORKDIR = joinpath(
            splitdir( dirname( FILEBRW ) )[ 1 ], splitext( basename( FILEBRW ) )[ 1 ] );
    end
    Result = [ ];
    for ( root, dirs, files ) in walkdir( WORKDIR )
        for file in files
            push!( Result, joinpath( root, file ) ); # path to files
        end
    end
    Result = Result[ occursin.( thing, Result ) ];
    if length( Result ) == 1
        Result = Result[ 1 ];
    end
    return Result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    RemezBPF( Variables::Dict, channel::Vector{Float64}, lF = 300,
        HF = 3000, fac = 10 )
        using DSP
"""
function RemezBPF( Variables::Dict, channel::Vector{Float64}, lF = 300,
        HF = 3000, fac = 10 )
    SamplingRate = Variables[ "SamplingRate" ];
    NYQ = floor( Int, SamplingRate / 2 );
    order = Int( floor( ( SamplingRate / lF ) / 5 ) );
    bpass = remez(
        ( order + 1 ), [ ( 0, lF - fac ) => 0, ( lF, HF ) => 1, ( HF + fac, NYQ ) => 0 ],
            Hz = SamplingRate );
    Filtered = filtfilt( DSP.Filters.PolynomialRatio( bpass, [ 1.0 ] ), channel );
    return Filtered
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    RemezBSF( Variables::Dict, channel::Vector{Float64}, lF = 60,
        HF = 60, fac = 5 )
        using DSP
"""
function RemezBSF( Variables::Dict, channel::Vector{Float64}, lF = 60,
        HF = 60, fac = 5 )
    SamplingRate = Variables[ "SamplingRate" ];
    NYQ = floor( Int, SamplingRate / 2 );
    order = Int( floor( ( SamplingRate / lF ) / 5 ) );
    bstop = remez(
            ( order + 1 ), [ ( 0, lF - fac ) => 1, ( lF, HF ) => 0, ( HF + fac, NYQ ) => 1 ],
                Hz = SamplingRate );
    Filtered = filtfilt( DSP.Filters.PolynomialRatio( bstop, [ 1.0 ] ), channel );
    return Filtered
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    MatrixFilter( Variables::Dict, data::Matrix{Float64}, lF = 300,
        HF = 3000, fac = 10 ) -> datafilt::Matrix{Float64}
        using RemezBPF, RemezBSF
"""
function MatrixFilter( Variables::Dict, data::Matrix{Float64}, lF,
        HF, fac, method::String )
    datafilt = copy( data );
    if method == "BPF"
        for k = 1:size( data, 1 )
            channel = data[ k, : ];
            Filtered = RemezBPF( Variables, channel, lF, HF, fac );
            datafilt[ k, : ] = Filtered;
        end
    elseif method == "BSF"
        for k = 1:size( data, 1 )
            channel = data[ k, : ];
            Filtered = RemezBSF( Variables, channel, lF, HF, fac );
            datafilt[ k, : ] = Filtered;
        end
    end
    return datafilt
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
Donoho( x ) = ( median( abs.( x ) ) / 0.6745 );
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
Energy( s, SamplingRate ) = sum( abs2, s ) / SamplingRate;
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
SigmaData( data ) = sqrt( 2*log( length( data ) ) );
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function OutlayerME( var1 )
    m = median( var1 );
    me, er = jackknife( identity, var1 );
    SePasa = vcat( findall( var1 .> ( m + er ) ), findall( var1 .< ( m - er ) ) );
    return SePasa
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function FindAllThrs( aux )
    a = zeros( 10 );
    a[ 1 ] = HistogramThresholding.find_threshold( aux, Balanced( ) );
    a[ 2 ] = HistogramThresholding.find_threshold( aux, Entropy( ) );
    a[ 3 ] = HistogramThresholding.find_threshold( aux, Intermodes( ) );
    a[ 4 ] = HistogramThresholding.find_threshold( aux, MinimumError( ) );
    a[ 5 ] = HistogramThresholding.find_threshold( aux, MinimumIntermodes( ) );
    a[ 6 ] = HistogramThresholding.find_threshold( aux, Moments( ) );
    a[ 7 ] = HistogramThresholding.find_threshold( aux, Otsu( ) );
    a[ 8 ] = HistogramThresholding.find_threshold( aux, UnimodalRosin( ) );
    a[ 9 ] = HistogramThresholding.find_threshold( aux, Yen( ) );
    a[ 10 ] = Donoho( aux )*SigmaData( aux );
    return a
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    LoadDict( s::String ) -> dicty::Dict
        using JLD
"""
function LoadDict( s::String )
    dicty = load( s );
    x = keys( dicty );
    if length( x ) == 1
        OK = collect( x )[ 1 ];
        dicty = dicty[ OK ];
    end
    return dicty
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function AllMs( Metrics, data, Parameters, Empties = [ ] );
    data = round.( Float64.( data ), digits = 3 );
    R = Array{ Real }( undef, size( data, 1 ), length( Metrics ) ); fill!( R, 0 );
    channels = collect( 1 : size( data, 1 ) );
    constant = sqrt( 2 * log( size( data, 2 ) ) );
    SpikeIndexes = [ ];
    for m in Metrics
        if m == 1; R[ :, m ] = UniqueCount( data ); end
        if m == 2; R[ :, m ] = std( data, dims = 2 ); end
        if m == 3; R[ :, m ] = stdΔV( Parameters, data, Parameters[ "ΔT" ] ); end
        if m == 4
            [ R[ k, m ] = Donoho( data[ k, : ] ) * constant for k in channels ];
        end
        if m == 5
            for k in channels
                _, R[ k, m ] = jackknife( identity, data[ k, : ] );
            end
        end
        if m == 6;
            IndexesSA = GetIndexes( Parameters, data, "Static", "All", Empties );
            R[ :, m ] = length.( IndexesSA );
            push!( SpikeIndexes, IndexesSA );
        end
        if m == 7;
            IndexesDW = GetIndexes( Parameters, data, "Dynamic", "Windows", Empties );
            R[ :, m ] = length.( IndexesDW );
            push!( SpikeIndexes, IndexesDW );
        end
        if m == 8; R[ :, m ] = median( data, dims = 2 ); end
        if m == 9; R[ :, m ] = mean( data, dims = 2 ); end
        if m == 10
            [ R[ k, m ] = sum( abs.( zscore( data[ k, : ] ) ) ) for k in channels ];
        end
    end
    return R, SpikeIndexes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function GetIndexes( Parameters,
                    data, ThrType = "Dynamic", Detection = "Windows", Empties = [ ] )
    if ThrType == "Static"
        SpikeIndexes = SpikeIndexesSA( data, Empties, Parameters );
    elseif ThrType == "Dynamic"
        if Detection == "Windows"
            SpikeIndexes = SpikeIndexesDW( data, Empties, Parameters );
        elseif Detection == "All";
            SpikeIndexes = SpikeIndexesDA( data, Empties, Parameters );
        end
    end
    return SpikeIndexes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function SpikeIndexesSA( data, Empties, Parameters )
    nChs, nFrs = size( data );
    ValidChannels = setdiff( 1:nChs, Empties );
    stokes = ms2Frs( Parameters, Parameters[ "stokes" ] );
    thr = Parameters[ "thr" ];
    SpikeIndexes = Array{ Any }( undef, nChs ); fill!( SpikeIndexes, [ ] );
    for vch in ValidChannels;
        signal = data[ vch, : ];
        I = findall( signal .<= thr );
        sobres = 1
        while !isempty( sobres )
            sobres = findall( diff( I ) .<= stokes );
            pares = hcat( I[ sobres ], I[ sobres .+ 1 ] );
            _, B = findmax( signal[ pares ], dims = 2 );
            I = I[ I .∉ [ pares[ B ] ] ];
        end
        SpikeIndexes[ vch ] = I;
    end
    return SpikeIndexes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function SpikeIndexesDW( data, Empties, Parameters )
    nChs, nFrs = size( data );
    ValidChannels = setdiff( 1:nChs, Empties );
    bit = ms2Frs( Parameters, Parameters[ "bit" ] );
    window = ms2Frs( Parameters, Parameters[ "window" ] );
    stokes = ms2Frs( Parameters, Parameters[ "stokes" ] )
    SpikeIndexes = Array{ Any }( undef, nChs ); fill!( SpikeIndexes, [ ] )
    for vch in ValidChannels;
        signal = data[ vch, : ];
        Is = collect( 1:bit:length( signal ) );
        aux = [ ];
        for I in Is
            E = I + window - 1;
            if E > length( signal )
                E = length( signal );
            end
            segment = signal[ I : E ];
            thr = -1 * Donoho( segment ) * SigmaData( segment );
            push!( aux, findall( segment .<= thr ) .+ ( I - 1 ) );
        end
        I = sort( unique( vcat( aux... ) ) );
        sobres = 1;
        while !isempty( sobres )
            sobres = findall( diff( I ) .<= stokes );
            pares = hcat( I[ sobres ], I[ sobres .+ 1 ] );
            _, B = findmax( signal[ pares ], dims = 2 );
            I = I[ I .∉ [ pares[ B ] ] ];
        end
        SpikeIndexes[ vch ] = I;
    end
    return SpikeIndexes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function SpikeIndexesDA( data, Empties, Parameters )
    nChs, nFrs = size( data );
    ValidChannels = setdiff( 1:nChs, Empties );
    stokes = ms2Frs( Parameters, Parameters[ "stokes" ] );
    SpikeIndexes = Array{ Any }( undef, nChs ); fill!( SpikeIndexes, [ ] )
    for vch in ValidChannels;
        signal = data[ vch, : ];
        thr = -1 * Donoho( signal ) * SigmaData( signal );
        I = findall( signal .<= thr );
        sobres = 1
        while !isempty( sobres )
            sobres = findall( diff( I ) .<= stokes );
            pares = hcat( I[ sobres ], I[ sobres .+ 1 ] );
            _, B = findmax( signal[ pares ], dims = 2 );
            I = I[ I .∉ [ pares[ B ] ] ];
        end
        SpikeIndexes[ vch ] = I;
    end
    return SpikeIndexes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    PFT_mt(  Variables::Dict, data::Matrix{Float64} ) -> S::Vector{Any}
        using SpectroMT, Spect2Ss
"""
function PFT_mt( Variables::Dict, data::Matrix{Float64} )
    Spect = SpectroMT( Variables, data );
    S = Spect2Ss( Spect );
    return S
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    Spect2Ss( Spect ) -> Ss::Vector{Any}
        using DSP
"""
function Spect2Ss( Spect )
    Ss = [ ];
    for sp in Spect
        T = collect( time( sp ) ); P = power( sp ); F = freq( sp );
        F = vcat( NaN, F );
        S = hcat( F, vcat( reshape( T, 1, length( T ) ), P ) );
        S = round.( S, digits = 4 ); S = Float16.( S );
        push!( Ss, S )
    end
    return Ss
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    SpectroMT( Variables::Dict, data::Matrix{Float64} )
        -> R::Array{DSP.Periodograms.Spectrogram}
        using DSP
"""
function SpectroMT( Variables::Dict, data::Matrix{Float64} )
    SamplingRate = Variables[ "SamplingRate" ];
    nChs, nFrs = size( data );
    channels = collect( 1:nChs );
    R = Array{DSP.Periodograms.Spectrogram}( undef, nChs );
    for ch in channels
        R[ ch ] = DSP.mt_spectrogram( data[ ch, : ], fs = floor( Int, SamplingRate ) );
    end
    return R
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function UnNormGaussK( x, sigma )
    return exp( ( -x * x ) / ( 2 * sigma ) );
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function UnNormGaussI( x, sigma )
    cte = ( 1 / ( sigma * sqrt( 2 * pi ) ) );
    return cte * exp( -( ( x ^ 2 ) / ( 2 * ( sigma ^ 2 ) ) ) );
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function GaussSmoothTemporal( Datos::Array, Sigma = 3, Gauss = "I" )
    # Un suavizado Gaussiano temporal.
    # Esto es escencialmente un filtro pasabajos.
    # Depende implicitamente de la frecuencia de muestreo.
    # sigma esta medido en pixeles, es la desviacion estandar de nuestro kernel.
    # El A de nuestra ventana seran 3 * sigma
    A = ceil( Sigma * 3 );
    B = ones( A );
    result = zeros( size( Datos ) );
    datosB = vcat( B * Datos[ 1 ], Datos, B * Datos[ end ] );
    if Gauss == "K"
        kernel = map( x -> UnNormGaussK( x, Sigma ), collect( -A : A ) );
        kernel = kernel / ( sum( kernel ) );
    elseif Gauss == "I"
        kernel = map( x -> UnNormGaussI( x, Sigma ), collect( -A : A ) );
        kernel = kernel / ( sum( kernel ) );
    end
    for t = ( A + 1 ) : ( length( Datos ) + A )
        result[ t - A ] = sum(
            datosB[ ( t - A ) : ( t + A ) ] .* kernel );
    end
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function DiscreteLaplacian( Data::Array )
    ( J, K ) = size( Data );
    I = reshape( Data[ 1, : ], ( 1, K ) );
    R = reshape( Data[ end, : ], ( 1, K ) );
    # Padding con copia de los datos
    Data = vcat( I, Data, R );
    Data = hcat( Data[ :, 1 ], Data, Data[ :, end ] );
    J, K = size( Data );
    Fx = Array{ Float64 }( undef, 3, 3 );
    result = zeros( size( Data ) );
    # calcular el CSD aplicando el Kernel laplaciano a cada celda más su 8-Vecindad y suma
    # todos los resultados como valor de la celda
    for j = ( 2 : ( J - 1 ) ), k = ( 2 : ( K - 1 ) )
        Fx = Data[ ( j - 1 ) : ( j + 1 ), ( k - 1 ) : ( k + 1 ) ];
        result[ j, k ] = sum( LaplacianKernel .* Fx );
    end
    # Crop the borders
    result = result[ 2 : ( end - 1 ) , 2 : ( end - 1 ) ]
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function GaussianBlur( Data::Array )
    result = zeros( size( Data ) );
    ( J, K ) = size( Data );
    # Padding con copia de los datos
    A = reshape( Data[ 1, : ], ( 1, K ) );
    B = reshape( Data[ end, : ], ( 1, K ) );
    C = vcat( A, A, A );
    D = vcat( B, B, B );
    Data = vcat( C, Data, D );
    for j = 1 : 3
        Data = hcat( Data[ :, 1 ], Data, Data[ :, end ] );
    end
    for j = ( 4 : ( J + 3 ) ), k = ( 4 : ( K + 3 ) )
        bin = Data[ ( j - 3 ) : ( j + 3 ), ( k - 3 ) : ( k + 3 ) ];
        result[ ( j - 3 ), ( k - 3 ) ] = sum( GaussianKernel .* bin );
    end
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function CSD( Data::Array, sigma = 3, Gauss = "I" )
    nChs, nFrs = size( Data );
    side = Int( sqrt( nChs ) );
    Data = reshape( Data, side, side, nFrs );
    ( J, K, L )  = size( Data );
    # We apply a Temporal Gaussian smoothing ( this greatly affects the animations )
    BINSMOOTH = zeros( J, K, L );
    for j = 1 : J, k = 1 : K
        BINSMOOTH[ j, k, : ] = GaussSmoothTemporal( vec( Data[ j, k, : ] ), sigma, Gauss );
    end
    ∇ = zeros( J, K, L );
    # We spatially smooth the LFP with a two-dimensional Gaussian filter.
    # Later we obtain the dCSD.
    for l = 1 : L
        ∇[ :, :, l ] = DiscreteLaplacian( GaussianBlur( BINSMOOTH[ :, :, l ] ) );
    end
    ∇ = -1 .* ∇;
    return ∇
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    CMs( csd::Array, t0::Int, tN::Int, ϵ::Real, minChans::Int ) -> CMN, CMP
"""
function CMs( csd::Array, t0::Int, tN::Int, ϵ::Real, minChans::Int )
    ( A, B, _ ) = size( csd );
    CMP = Dict{ Int, Array }( );
    CMN = Dict{ Int, Array }( );
    for t = t0:tN
        NegChans = Array{ Int16 }[ ];
        PosChans = Array{ Int16 }[ ];
        SpikeCountPositivo = zeros( A, B );
        SpikeCountNegativo = zeros( A, B );
        for a = ( 1 : A ), b = ( 1 : B )
            if( csd[ a, b, t ] < ( -1 * ϵ ) )
                push!( NegChans, [ a, b ] );
                SpikeCountNegativo[ a, b ] += 1;
            elseif( csd[ a, b, t ] > ϵ )
                push!( PosChans, [ a, b ] );
                SpikeCountPositivo[ a, b ] += 1;
            end
        end
        CMN[ t ] = CentrosMasa( csd, t, NegChans, minChans );
        CMP[ t ] = CentrosMasa( csd, t, PosChans, minChans );
    end
    return CMN, CMP
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
"""
    CentrosMasa( CSDA, canales, minchannels = 3 )
            CM = [ x y masa ];
"""
function CentrosMasa( CSDA::Array, t::Real, canales::Vector, minchannels = 3 )
    grupos = ComponentesSP( canales );
    centros = [ 0 0 0 ];
    for p in grupos
        mu = length( p );
        if mu >= minchannels
            masa = 0.00;
            x = 0.00;
            y = 0.00;
            for q in p
                masalocal = CSDA[ q[ 1 ], q[ 2 ], t ];
                masa += masalocal;
                x += q[ 2 ] * masalocal;
                y += q[ 1 ] * masalocal;
            end
            x/=masa
            y/=masa
            CM = [ x y masa ];
            centros = vcat( centros, CM )
        end
    end
    centros = centros[ ( 2 : end ), : ];
    return centros
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function TiraOrillas( Puntos::Set )
    # Descarta lo que se sale de la malla de electrodos
    result = Set( [ ] )
    for p in Puntos
        if !( p[ 1 ] == 0 || p[ 2 ] == 0 || p[ 1 ] == 65 || p[ 2 ] == 65 )
            push!( result, p );
        end
    end
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function vecindad8( punto::Array )
    # La ocho-vecindad de un punto en una malla cuadrada.
    j = punto[ 1 ];
    k = punto[ 2 ];
    result = Set{ Array{ Int64, 1 } }( );
    push!( result, [ j - 1, k - 1 ] );
    push!( result, [ j - 1, k ] );
    push!( result, [ j - 1, k + 1 ] );
    push!( result, [ j, k - 1 ] );
    push!( result, [ j, k + 1 ] );
    push!( result, [ j + 1, k - 1 ] );
    push!( result, [ j + 1, k ] );
    push!( result, [ j + 1, k + 1 ] );
    result = TiraOrillas( result )
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function ComponentesSP( DatosSignados::Array )
    #Single pass method for Disjoint Components.
    lista = copy( DatosSignados );
    componentes = Set{ Any }( );
    while ( length( lista ) != 0 )
        x = pop!( lista ); #arranca el ULTIMO elemento de la lista
        listaprofundeza = Array{ Int64 }[ ];
        componentecurlab = Array{ Int64 }[ ];
        push!( listaprofundeza, x ); #Pone elementos al FINAL de la lista
        push!( componentecurlab, x );
        profundidad = 0;
        while ( ( length( listaprofundeza ) != 0 ) && profundidad < 1000 )
            y = pop!( listaprofundeza );
            for v in vecindad8( y )
                if in( v, lista ) # A: Si v está en la lista
                    deleteat!( lista, indexin( Any[ v ], lista ) );
                    push!( listaprofundeza, v );
                    profundidad += 1;
                    push!( componentecurlab, v );
                end
            end
        end
        push!( componentes, componentecurlab );
    end
    return componentes
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #

end # module
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
