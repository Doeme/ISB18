<map version="1.0.1">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1523621441590" ID="ID_46461404" MODIFIED="1523621449727" TEXT="final report">
<node CREATED="1523621472908" ID="ID_1254202059" MODIFIED="1523621477101" POSITION="right" TEXT="Introduction">
<node CREATED="1523622188024" ID="ID_1112896377" MODIFIED="1523622193091" TEXT="basic problem:">
<node CREATED="1523622195172" ID="ID_190233608" MODIFIED="1523622232785" TEXT="in alcohole fermentation in food industry (wine, beer)"/>
<node CREATED="1523625370979" ID="ID_1714801516" MODIFIED="1523625376778" TEXT="take beer brewing as an example"/>
<node CREATED="1523622236244" ID="ID_919210626" MODIFIED="1523622252152" TEXT="high quality requirements"/>
<node CREATED="1523622252453" ID="ID_919038601" MODIFIED="1523622271127" TEXT="many restrictions how to influence the fermentation"/>
<node CREATED="1523622278606" ID="ID_101013873" MODIFIED="1523622289003" TEXT="leads to a complex optimization problem"/>
<node CREATED="1523622296875" ID="ID_1577525145" MODIFIED="1523622305155" TEXT="high effort in laboratory"/>
</node>
<node CREATED="1523622316797" ID="ID_324441527" MODIFIED="1523622320234" TEXT="basic solution">
<node CREATED="1523624072474" ID="ID_217489868" MODIFIED="1523624083773" TEXT="simulation of the fermentation process"/>
</node>
<node CREATED="1523624098854" ID="ID_181795023" MODIFIED="1523624153775" TEXT="Give an overview over existing approaches (general simulation approaches)"/>
<node CREATED="1523627809096" ID="ID_1586135011" MODIFIED="1523627825059" TEXT="why DFBA?"/>
<node CREATED="1523621596699" ID="ID_1724314888" MODIFIED="1523624944673" TEXT="Using GEMs and DFBA in alcohole fermentation is a modern tool to control the process">
<node CREATED="1523622175744" ID="ID_502350890" MODIFIED="1523622178189" TEXT="example wine"/>
</node>
<node CREATED="1523624177568" ID="ID_355010992" MODIFIED="1523624935586" TEXT="Give an overview over interesting implementations"/>
<node CREATED="1523625343149" ID="ID_1939763195" MODIFIED="1523625396225" TEXT="discussion: which implementation fits best to this problem"/>
<node CREATED="1523621602700" ID="ID_1025772164" MODIFIED="1523625872814" TEXT="Which results does it produce?"/>
<node CREATED="1523621607137" ID="ID_157481879" MODIFIED="1523621786122" TEXT="How does it do that?"/>
</node>
<node CREATED="1523621477528" ID="ID_1383421190" MODIFIED="1523621487179" POSITION="right" TEXT="Methods">
<node CREATED="1523626581936" ID="ID_688850775" MODIFIED="1523627045345" TEXT="Genome-Scale Models (GEMs)"/>
<node CREATED="1523626453603" ID="ID_23882357" MODIFIED="1523626470592" TEXT="Flux Balance Analysis (FBA)"/>
<node CREATED="1523626611806" ID="ID_1595821955" MODIFIED="1523627016567">
<richcontent TYPE="NODE"><html>
  <head>
    
  </head>
  <body>
    <p>
      Dynamic Flux Balance Analysis (DFBA)
    </p>
    <p>
      (our implementation)
    </p>
  </body>
</html>
</richcontent>
<node CREATED="1523627114019" ID="ID_835940491" MODIFIED="1523627128381" TEXT="our implementation orientates on DMMM"/>
<node CREATED="1523627128894" ID="ID_738096624" MODIFIED="1523627140102" TEXT="list up some differences to DMMM">
<node CREATED="1523627142560" ID="ID_1406425672" MODIFIED="1523627155972" TEXT="(see table in report #2)"/>
</node>
<node CREATED="1523627394177" ID="ID_1900593750" MODIFIED="1523627418262" TEXT="short description of the algorithm"/>
</node>
<node CREATED="1523627432228" ID="ID_827184541" MODIFIED="1523627440078" TEXT="Bottlenecks"/>
<node CREATED="1523627352961" ID="ID_1500317686" MODIFIED="1523627362944" TEXT="Simulation Setup">
<node CREATED="1523627365600" ID="ID_1750906473" MODIFIED="1523627369495" TEXT="used GEMs"/>
<node CREATED="1523627370925" ID="ID_405200565" MODIFIED="1523627380852" TEXT="constants"/>
</node>
</node>
<node CREATED="1523621487455" ID="ID_14969992" MODIFIED="1523621491036" POSITION="right" TEXT="Results"/>
<node CREATED="1523621491367" ID="ID_1087729205" MODIFIED="1523621493468" POSITION="right" TEXT="Conclusions">
<node CREATED="1523627451390" ID="ID_1656285519" MODIFIED="1523627708217" TEXT="How can the simulation be enhanced?">
<node CREATED="1523627457760" ID="ID_1493169712" MODIFIED="1523627510864" TEXT="how to get the right infput/output fluxes for each time point? (FBA/FVA)"/>
<node CREATED="1523627513012" ID="ID_536727355" MODIFIED="1523627517609" TEXT="how to model the mortality"/>
<node CREATED="1523627574553" ID="ID_1847290381" MODIFIED="1523627617225" TEXT="trade between many data points (accuracy) and calculation effort (time consumption)"/>
<node CREATED="1523627517875" ID="ID_1556904958" MODIFIED="1523628147299" TEXT="no toxicity"/>
<node CREATED="1523627530071" ID="ID_1887405089" MODIFIED="1523627549423" TEXT="no temperatur"/>
</node>
<node CREATED="1523628172435" ID="ID_1178464299" MODIFIED="1523628175385" TEXT="Outlook">
<node CREATED="1523628182078" ID="ID_924714734" MODIFIED="1523628186770" TEXT="next steps?">
<node CREATED="1523628203345" ID="ID_436939676" MODIFIED="1523628238503" TEXT="metabolite austausch mit umgebung:">
<node CREATED="1523628240691" ID="ID_1926795350" MODIFIED="1523628260988" TEXT="input flux of oxygen via solution surface"/>
</node>
<node CREATED="1523628264914" ID="ID_1984953785" MODIFIED="1523628358184" TEXT="verwertung der bestandteile von abgestorbenen bakterien in der L&#xf6;sung"/>
</node>
</node>
</node>
</node>
</map>
