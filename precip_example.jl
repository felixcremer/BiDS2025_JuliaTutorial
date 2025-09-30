using YAXArrays, NetCDF
import DimensionalData as DD
import Dates
using CairoMakie

#Open the time zones
tzfile = "/Net/Groups/BGI/work_2/scratch/mweynants/ARCEME/timezone_offset.nc"
tz = open_dataset(tzfile)
tz.offset.data |> unique
heatmap(tz.offset,colorrange=(-12,12))

#Open the precip data
p_era = "/Net/Groups/data_BGC/era5_land/e1/0d10_hourly/tp"
allmonths = [filter(endswith(".nc"),sort(readdir(joinpath(p_era,"$yr"),join=true))) for yr in 1990:2023]
allmonths = reduce(vcat,allmonths)
_precip = open_mfdataset(DD.DimArray(allmonths,YAXArrays.time()),skip_keys=(:number,)).tp
precip = cat(_precip,_precip,dims=longitude(-360.0:0.1:359.9))[longitude=-180..179.9]

#Plot a single time step
heatmap(precip[time=Near(DateTime(2015,6,23,13,0,0))].^0.25)

#Look at a time series
subset = longitude(Near(11.0)), latitude(Near(51.0)), YAXArrays.time(DateTime(2015,6,22,1)..DateTime(2015,6,25))

data_subset = precip[subset...][:]
lines(data_subset)

#Compute instantaneous precip 
struct DailyDiffVector{T,P<:AbstractVector{T}} <: AbstractVector{T}
    parent::P
end
Base.size(a::DailyDiffVector) = (length(a.parent)-1,)
function Base.getindex(a::DailyDiffVector,i::Int) 
    if mod(i-1,24) == 0 
        # We are in a first hour of day
        a.parent[i]
    else
        # We need to take diff to previous
        (a.parent[i]-a.parent[i-1])
    end
end

precip_inst = DailyDiffVector(data_subset.data)
lines(precip_inst)

#Applying the time zone to get the daily sum for a calendar day
function precipsum(xin,tzoffset)
    
    if ismissing(xin) || tzoffset < -20 || tzoffset > 20 
        return NaN
    end
    aoi = (25:48) .- round(Int,tzoffset)
    b = DailyDiffVector(xin)
    current_day = view(b,aoi)
    dailyprecip = sum(current_day)
    return Float64(dailyprecip)
end

w = MovingIntervals(:closed, :open, left=DateTime(1990,1,1,1):Dates.Day(1):DateTime(2023,12,29,1),width=Dates.Day(3))

precipwindows = windows(precip,:time=>w)

xprecipstats = XFunction(precipsum,inplace=false)
dailyprecip = xprecipstats.(precipwindows,tz.offset)

firstday = dailyprecip[longitude=-5.0..35.0,latitude=35.0..65.0,time=1]

firstday_computed = firstday[:,:]

heatmap(firstday_computed.^(0.25))


outds = Dataset(prec = dailyprecip)

using Distributed, SlurmClusterManager

addprocs(SlurmClusterManager())
@everywhere using YAXArrays, NetCDF
@everywhere include("precipstats.jl")

compute_to_zarr(outds,"./dailyprecip.zarr",overwrite=true)

rmprocs(workers())
