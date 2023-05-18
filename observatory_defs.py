from Observatory import Observatory
import Utilities

from Swope import Swope
from Nickel import Nickel
from Thacher import Thacher
from T80S import T80S
from Blanco import Blanco

lco = Observatory(
        name="LCO",
        lon="-70.6915",
        lat="-29.0182",
        elevation=2402,
        horizon="-12",
        telescopes={"Swope":Swope()},
        utc_offset=Utilities.lco_clt_utc_offset,
        utc_offset_name="CLST",
)

lick = Observatory(
        name="Lick",
        lon="-121.6429",
        lat="37.3414",
        elevation=1283,
        horizon="-12",
        telescopes={"Nickel":Nickel()},
        utc_offset=Utilities.lick_pst_utc_offset,
        utc_offset_name="PST",
)

thacher = Observatory(
        name="Thacher",
        lon="-121.6431",
        lat="34.46479",
        elevation=630.0,
        horizon="-12",
        telescopes={"Thacher":Thacher()},
        utc_offset=Utilities.lick_pst_utc_offset,
        utc_offset_name="PST",
)

ctio = Observatory(
        name="CTIO",
        lon="-70.8035",
        lat="-30.1732",
        #lon="-70.6915",
        #lat="-29.0182",
        elevation=2207.0,
        horizon="-12",
        telescopes={"T80S":T80S(),"Blanco":Blanco()},
        utc_offset=Utilities.lco_clt_utc_offset,
        utc_offset_name="CLST",
)

observatories = {"LCO":lco, "Lick":lick, "Thacher":thacher,"CTIO":ctio}
