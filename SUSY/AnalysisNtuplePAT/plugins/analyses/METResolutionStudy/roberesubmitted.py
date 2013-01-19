#/usr/bin/python

import sys,os
import string
import re

killed = [
    "34736.0",
    "34738.0",
    "34739.0",
    "34740.0",
    "34741.0",
    "34742.0",
    "34743.0",
    "34745.0",
    "34746.0",
    "34747.0",
    "34749.0",
    "34750.0",
    "34751.0",
    "34753.0",
    "34754.0",
    "34755.0",
    "34756.0",
    "34757.0",
    "34758.0",
    "34759.0",
    "34760.0",
    "34762.0",
    "34763.0",
    "34764.0",
    "34765.0",
    "34766.0",
    "34767.0",
    "34768.0",
    "34769.0",
    "34770.0",
    "34771.0",
    "34772.0",
    "34773.0",
    "34774.0",
    "34775.0",
    "34776.0",
    "34777.0",
    "34778.0",
    "34779.0",
    "34780.0",
    "34781.0",
    "34782.0",
    "34785.0",
    "34786.0",
    "34787.0",
    "34788.0",
    "34789.0",
    "34790.0",
    "34791.0",
    "34792.0",
    "34793.0",
    "34795.0",
    "34796.0",
    "34797.0",
    "34798.0",
    "34799.0",
    "34800.0",
    "34801.0",
    "34802.0",
    "34803.0",
    "34804.0",
    "34805.0",
    "34806.0",
    "34807.0",
    "34808.0",
    "34809.0",
    "34810.0",
    "34811.0",
    "34812.0",
    "34813.0",
    "34816.0",
    "34817.0",
    "34818.0",
    "34819.0",
    "34820.0",
    "34821.0",
    "34822.0",
    "34823.0",
    "34824.0",
    "34826.0",
    "34827.0",
    "34828.0",
    "34829.0",
    "34830.0",
    "34831.0",
    "34832.0",
    "34833.0",
    "34834.0",
    "34836.0",
    "34837.0",
    "34838.0",
    "34839.0",
    "34840.0",
    "34842.0",
    "34843.0",
    "34844.0",
    "34845.0",
    "34846.0",
    "34848.0",
    "34849.0",
    "34850.0",
    "34851.0",
    "34852.0",
    "34855.0",
    "34857.0",
    "34858.0",
    "34859.0",
    "34860.0",
    "34861.0",
    "34864.0",
    "34865.0",
    "34866.0",
    "34867.0",
    "34868.0",
    "34869.0",
    "34870.0",
    "34871.0",
    "34872.0",
    "34873.0",
    "34874.0",
    "34875.0",
    "34876.0",
    "34877.0",
    "34878.0",
    "34879.0",
    "34880.0",
    "34881.0",
    "34882.0",
    "34883.0",
    "34884.0",
    "34885.0",
    "34886.0",
    "34887.0",
    "34889.0",
    "34891.0",
    "34892.0",
    "34893.0",
    "34894.0",
    "34895.0",
    "34896.0",
    "34898.0",
    "34899.0",
    "34900.0",
    "34902.0",
    "34903.0",
    "34904.0",
    "34906.0",
    "34907.0",
    "34910.0",
    "34911.0",
    "34912.0",
    "34913.0",
    "34916.0",
    "34919.0",
    "34920.0",
    "34921.0",
    "34922.0",
    "34923.0",
    "34924.0",
    "34925.0",
    "34926.0",
    "34927.0",
    "34928.0",
    "34929.0",
    "34930.0",
    "34932.0",
    "34933.0",
    "34934.0",
    "34936.0",
    "34937.0",
    "34939.0",
    "34942.0",
    "34943.0",
    "34946.0",
    "34947.0",
    "34948.0",
    "34949.0",
    "34950.0",
    "34951.0",
    "34952.0",
    "34953.0",
    "34954.0",
    "34955.0",
    "34956.0",
    "34957.0",
    "34958.0",
    "34959.0",
    "34960.0",
    "34961.0",
    "34963.0",
    "34968.0",
    "34969.0",
    "34970.0",
    "34972.0",
    "34973.0",
    "34975.0",
    "34976.0",
    "34977.0",
    "34978.0",
    "34979.0",
    "34980.0",
    "34982.0",
    "34983.0",
    "34985.0",
    "34986.0",
    "34987.0",
    "34989.0",
    "34990.0",
    "34991.0",
    "34992.0",
    "34993.0",
    "34994.0",
    "34995.0",
    "34996.0",
    "34997.0",
    "34998.0",
    "34999.0",
    "35000.0",
    "35001.0",
    "35002.0",
    "35003.0",
    "35005.0",
    "35007.0",
    "35008.0",
    "35009.0",
    "35010.0",
    "35011.0",
    "35014.0",
    "35015.0",
    "35016.0",
    "35017.0",
    "35018.0",
    "35019.0",
    "35020.0",
    "35022.0",
    "35023.0",
    "35024.0",
    "35025.0",
    "35026.0",
    "35027.0",
    "35028.0",
    "35029.0",
    "35030.0",
    "35031.0",
    "35034.0",
    "35035.0",
    "35036.0",
    "35037.0",
    "35038.0",
    "35039.0",
    "35040.0",
    "35041.0",
    "35042.0",
    "35043.0",
    "35044.0",
    "35045.0",
    "35046.0",
    "35047.0",
    "35049.0",
    "35050.0",
    "35051.0",
    "35052.0",
    "35053.0",
    "35054.0",
    "35055.0",
    "35056.0",
    "35057.0",
    "35058.0",
    "35059.0",
    "35060.0",
    "35061.0",
    "35062.0",
    "35063.0",
    "35064.0",
    "35065.0",
    "35066.0",
    "35069.0",
    "35070.0",
    "35073.0",
    "35074.0",
    "35075.0",
    "35076.0",
    "35078.0",
    "35079.0",
    "35080.0",
    "35081.0",
    "35082.0",
    "35083.0",
    "35084.0",
    "35085.0",
    "35086.0",
    "35087.0",
    "35088.0",
    "35089.0",
    "35090.0",
    "35091.0",
    "35092.0",
    "35093.0",
    "35094.0",
    "35095.0",
    "35096.0",
    "35097.0",
    "35098.0",
    "35099.0",
    "35101.0",
    "35102.0",
    "35103.0",
    "35104.0",
    "35105.0",
    "35106.0",
    "35107.0",
    "35108.0",
    "35109.0",
    "35111.0",
    "35112.0",
    "35113.0",
    "35115.0",
    "35116.0",
    "35117.0",
    "35119.0",
    "35120.0",
    "35123.0",
    "35124.0",
    "35125.0",
    "35126.0",
    "35128.0",
    "35129.0",
    "35130.0",
    "35131.0",
    "35132.0",
    "35133.0",
    "35134.0",
    "35135.0",
    "35136.0",
    "35137.0",
    "35138.0",
    "35139.0",
    "35140.0",
    "35141.0",
    "35142.0",
    "35143.0",
    "35144.0",
    "35145.0",
    "35146.0",
    "35147.0",
    "35148.0",
    "35149.0",
    "35150.0",
    "35151.0",
    "35152.0",
    "35153.0",
    "35154.0",
    "35155.0",
    "35156.0",
    "35157.0",
    "35158.0",
    "35159.0",
    "35160.0",
    "35169.0",
    "35170.0",
    "35171.0",
    "35172.0",
    "35175.0",
    "35176.0",
    "35177.0",
    "35178.0",
    "35179.0",
    "35180.0",
    "35181.0",
    "35183.0",
    "35184.0",
    "35185.0",
    "35186.0",
    "35188.0",
    "35192.0",
    "35195.0",
    "35196.0",
    "35197.0",
    "35199.0",
    "35200.0",
    "35201.0",
    "35203.0",
    "35204.0",
    "35205.0",
    "35206.0",
    "35207.0",
    "35208.0",
    "35209.0",
    "35210.0",
    "35211.0",
    "35212.0",
    "35213.0",
    "35214.0",
    "35215.0",
    "35216.0",
    "35217.0",
    "35218.0",
    "35219.0",
    "35220.0",
    "35221.0",
    "35223.0",
    "35224.0",
    "35225.0",
    "35226.0",
    "35227.0",
    "35228.0",
    "35229.0",
    "35230.0",
    "35231.0",
    "35232.0",
    "35233.0",
    "35234.0",
    "35235.0",
    "35236.0",
    "35237.0",
    "35238.0",
    "35240.0",
    "35241.0",
    "35242.0",
    "35243.0",
    "35244.0",
    "35245.0",
    "35246.0",
    "35248.0",
    "35250.0",
    "35251.0",
    "35252.0",
    "35253.0",
    "35254.0",
    "35255.0",
    "35256.0",
    "35257.0",
    "35258.0",
    "35259.0",
    "35260.0",
    "35261.0",
    "35262.0",
    "35263.0",
    "35264.0",
    "35265.0",
    "35266.0",
    "35267.0",
    "35270.0",
    "35271.0",
    "35272.0",
    "35273.0",
    "35274.0",
    "35275.0",
    "35276.0",
    "35277.0",
    "35278.0",
    "35279.0",
    "35280.0",
    "35281.0",
    "35282.0",
    "35283.0",
    "35284.0",
    "35285.0",
    "35286.0",
    "35287.0",
    "35288.0",
    "35289.0",
    "35290.0",
    "35291.0",
    "35292.0",
    "35293.0",
    "35294.0",
    "35295.0",
    "35296.0",
    "35297.0",
    "35298.0",
    "35300.0",
    "35301.0",
    "35302.0",
    "35308.0",
    "35309.0",
    "35311.0",
    "35314.0",
    "35317.0",
    "35318.0",
    "35319.0",
    "35321.0",
    "35322.0",
    "35323.0",
    "35324.0",
    "35326.0",
    "35327.0",
    "35328.0",
    "35329.0",
    "35330.0",
    "35331.0",
    "35332.0",
    "35333.0",
    "35334.0",
    "35335.0",
    "35337.0",
    "35338.0",
    "35339.0",
    "35340.0",
    "35341.0",
    "35342.0",
    "35343.0",
    "35344.0",
    "35347.0",
    "35348.0",
    "35350.0",
    "35351.0",
    "35352.0",
    "35353.0",
    "35354.0",
    "35355.0",
    "35356.0",
    "35357.0",
    "35358.0",
    "35359.0",
    "35360.0",
    "35362.0",
    "35363.0",
    "35364.0",
    "35365.0",
    "35366.0",
    "35367.0",
    "35368.0",
    "35369.0",
    "35370.0",
    "35371.0",
    "35376.0",
    "35377.0",
    "35378.0",
    "35379.0",
    "35380.0",
    "35381.0",
    "35385.0",
    "35386.0",
    "35390.0",
    "35391.0",
    "35392.0",
    "35394.0",
    "35395.0",
    "35397.0",
    "35398.0",
    "35400.0",
    "35401.0",
    "35402.0",
    "35403.0",
    "35404.0",
    "35405.0",
    "35407.0",
    "35408.0",
    "35409.0",
    "35410.0",
    "35411.0",
    "35412.0",
    "35413.0",
    "35414.0",
    "35415.0",
    "35416.0",
    "35417.0",
    "35418.0",
    "35419.0",
    "35420.0",
    "35421.0",
    "35422.0",
    "35423.0",
    "35424.0",
    "35425.0",
    "35426.0",
    "35428.0",
    "35429.0",
    "35431.0",
    "35432.0",
    "35433.0",
    "35435.0",
    "35436.0",
    "35437.0",
    "35438.0",
    "35439.0",
    "35440.0",
    "35441.0",
    "35442.0",
    "35443.0",
    "35445.0",
    "35446.0",
    "35447.0",
    "35448.0",
    "35449.0",
    "35451.0",
    "35452.0",
    "35453.0",
    "35454.0",
    "35455.0",
    "35456.0",
    "35457.0",
    "35458.0",
    "35459.0",
    "35460.0",
    "35461.0",
    "35462.0",
    "35463.0",
    "35466.0",
    "35468.0",
    "35469.0",
    "35470.0",
    "35471.0",
    "35472.0",
    "35473.0",
    "35474.0",
    "35475.0",
    "35476.0",
    "35477.0",
    "35478.0",
    "35479.0",
    "35480.0",
    "35481.0",
    "35482.0",
    "35483.0",
    "35484.0",
    "35485.0",
    "35487.0",
    "35488.0",
    "35489.0",
    "35492.0",
    "35493.0",
    "35494.0",
    "35495.0",
    "35496.0",
    "35498.0",
    "35500.0",
    "35502.0",
    "35503.0"
]

for job in killed:
    cmd = "condor_rm "+job
    print(cmd)
    os.system(cmd)