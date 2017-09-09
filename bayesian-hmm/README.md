## A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging

- [A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging](http://homepages.inf.ed.ac.uk/sgwater/papers/acl07-bhmm.pdf)
- [実装について](http://musyoku.github.io/2017/01/28/A-Fully-Bayesian-Approach-to-Unsupervised-Part-of-Speech-Tagging/)

## 準備

### macOS

macOSの場合、PythonとBoostはともにbrewでインストールする必要があります。

#### Python 3のインストール

```
brew install python3
```

`PYTHONPATH`を変更する必要があるかもしれません。

#### Boostのインストール

```
brew install boost-python --with-python3
```

### Linux

#### Boostのインストール

```
./bootstrap.sh --with-libraries=python --with-python=python3 --with-python-root=YOUR_PYTHON_ROOT
./b2 python=3.6 -d2 -j4
```

### ビルド

以下のコマンドで`bhmm.so`が生成され、Pythonから利用できるようになります。

```
make install
```

`makefile`のBoostのincludeパスを自分の環境に合わせて書き換えてください。

### MeCabのインストール

```
pip install mecab-python3
```

## 学習

英語のテキストファイルの場合は以下のコマンドで学習できます。

```
python3 train_en.py  -f ../../text/ptb.txt -split 1 --start-temperature 1.5 --min-temperature 0.08 -tags 16 -alpha 1 -beta 1 -e 20000
```

日本語のテキストファイルの場合は以下のコマンドで学習できます。

```
python3 train_ja.py  -f ../../text/kokoro.txt -split 1 --start-temperature 1.5 --min-temperature 0.08 -tags 16 -alpha 1 -beta 1 -e 40000
```

## 結果の可視化

各予測タグとそれに割り当てられた単語を表示するには以下のコマンドを実行します。

```
python3 tags.py -n 200
```

ハイパーパラメータの値も表示されます。

混同行列をプロットするには以下のコマンドを実行します。

```
python3 plot_ja.py -f ../../text/kokoro.txt
```

学習時に使用したテキストファイルを指定します。

## 結果

### Penn Treebank

`ptb.txt`の学習結果です。

```
python train_en.py  -f ../../text/ptb.txt -split 1 --start-temperature 1.5 --min-temperature 0.08 -tags 16 -alpha 1 -beta 1 -e 20000
```

![result](https://raw.githubusercontent.com/musyoku/images/master/pos_tagging/bhmm/ptb.png)

```
[1]
say(9061) and(706) ago(403) think(353) earlier(344) add(331) note(269) believe(227) report(227) be(218) <unk>(197) know(192) inc(162) argue(129) suggest(123) show(121) estimate(120) Co(118) do(109) indicate(104) mean(101) find(100) Corp(83) tell(83) contend(83) feel(77) announce(76) predict(74) claim(74) insist(67) see(63) succeed(59) ask(58) agree(58) old(58) acknowledge(54) warn(51) clear(51) concede(50) whose(47) conclude(47) include(45) declare(45) ltd(44) complain(43) confirm(43) decide(42) call(41) name(41) concern(41) later(41) explain(40) fear(38) figure(37) post(37) disclose(36) realize(36) worry(36) assume(36) rule(35) hope(34) determine(34) recall(33) charge(32) where(32) admit(32) here(31) wonder(31) write(30) deny(30) corp(30) officer(29) allege(29) nor(27) maker(26) assert(26) or(25) question(25) maintain(25) expect(24) understand(24) doubt(24) plc(22) confident(22) president(21) judge(21) whom(21) suspect(21) hold(20) state(20) recommend(20) advise(19) before(18) sure(18) surprising(18) observe(18) discover(17) ensure(16) reply(16) calif(16) want(15) 
[2]
that(5148) but(3170) mr.(2378) Mr.(1623) and(1346) if(1199) when(923) because(619) as(607) so(395) while(393) what(358) even(325) which(309) however(288) though(287) now(257) although(255) where(212) then(200) whether(183) some(177) after(153) how(134) since(133) before(130) yesterday(124) meanwhile(122) mrs.(120) still(119) all(118) ms.(114) many(113) than(112) <unk>(108) president(101) yet(100) moreover(93) until(91) why(84) not(83) indeed(83) Ms.(81) also(75) once(75) unless(65) far(62) instead(60) like(57) maybe(56) first(52) today(52) such(52) judge(51) separately(51) thus(49) perhaps(46) do(45) qintex(45) Dr.(44) only(43) most(43) unlike(42) dr.(42) add(41) earlier(40) one(39) be(35) recently(35) make(35) currently(32) bond(30) Mrs.(30) both(28) about(28) well(26) later(26) whatever(24) nevertheless(24) gen.(24) sometimes(24) or(23) despite(23) already(23) just(22) note(22) consider(21) ask(21) philip(21) neither(21) finally(21) jack(21) unfortunately(20) time(19) fannie(19) can(17) bear(17) here(17) who(16) george(16) often(16) 
[3]
rise(889) price(719) net(675) stock(651) rate(611) sale(496) bond(475) fall(440) income(418) close(399) <unk>(395) revenue(389) profit(376) earning(344) average(325) increase(320) exchange(308) compare(295) share(278) gain(276) index(275) trading(268) total(267) third-quarter(261) loss(254) drop(249) interest(243) end(232) high(224) due(221) more(216) be(214) annual(210) yesterday(198) issue(191) composite(187) volume(183) yield(171) low(169) decline(165) value(165) estimate(153) operate(151) report(148) note(147) mortgage(131) jump(127) market(119) offering(112) offer(110) amount(109) security(107) cost(105) convertible(98) common(98) plunge(95) post(95) pretax(93) asset(89) late(87) growth(81) treasury(81) climb(81) charge(80) result(79) dividend(78) Friday(76) loan(74) tax(71) quarterly(70) fund(68) October(68) delivery(68) earn(66) because(66) debenture(66) grow(65) cash(65) base(63) 30-year(63) debt(62) September(62) dollar(61) slightly(61) surge(61) contract(60) future(60) Monday(59) subordinated(58) account(58) quote(57) return(55) as(55) face(54) preferred(54) unchanged(54) Nov.(54) range(53) coupon(53) trade(51) advance(50) 
[4]
<unk>(9437) and(1866) &(1388) .(1147) a(1033) president(854) corp.(817) co.(745) Inc.(706) chief(596) executive(591) chairman(527) vice(428) officer(337) inc.(315) an(309) mr.(285) john(282) general(271) Corp.(269) director(257) analyst(256) group(220) international(220) unit(209) first(202) Co.(199) security(192) james(188) motor(186) robert(185) financial(182) d(177) ##(163) former(163) capital(152) ltd.(152) david(149) the(135) senior(134) j(132) at(130) national(129) service(124) American(124) system(122) industry(121) new(119) William(119) lynch(116) manager(110) bank(109) t(106) computer(106) investment(106) associate(106) rep.(106) pacific(105) r(105) peter(103) richard(103) large(100) economist(98) l(95) electric(94) communication(94) insurance(92) management(90) boston(90) investor(89) Merrill(89) texas(88) operate(88) name(87) brother(86) de(84) market(81) attorney(80) say(78) research(77) trust(76) holding(76) e(73) plc(73) Salomon(72) g(72) smith(72) city(71) Co(70) b(70) poor(70) partner(69) standard(69) head(69) charles(69) Michael(68) chemical(68) goldman(67) c(67) subsidiary(66) lead(66) 
[5]
of(26152) in(18919) and(11417) for(9562) to(8335) on(5917) at(5003) by(4929) with(4900) as(4221) from(3720) be(2039) than(1608) about(1266) that(1153) or(1124) into(1014) after(990) over(878) <unk>(760) under(676) include(664) through(614) against(537) during(507) up(493) since(477) between(461) when(432) before(412) among(390) but(380) like(303) off(292) while(290) only(256) until(254) out(252) within(217) without(217) despite(198) reflect(176) around(171) via(170) follow(163) above(162) where(159) earlier(158) if(151) down(150) because(145) below(145) all(139) toward(130) too(120) call(100) give(97) last(94) across(90) involve(82) end(81) have(80) yesterday(79) near(79) behind(78) represent(72) so(68) beyond(66) just(63) take(62) say(62) outside(62) nearly(60) amid(60) throughout(58) rate(53) later(53) whose(51) indicate(51) back(51) less(49) hit(49) once(49) carry(48) corp.(46) cite(45) almost(45) along(44) leave(42) plus(41) use(39) how(39) time(38) much(37) exceed(37) fee(36) away(34) although(34) both(33) particularly(32) air(31) 
[6]
##(31302) $(8361) million(4982) share(2035) billion(2026) a(1837) about(933) cent(922) up(585) point(395) yen(364) down(356) oct.(355) day(304) yield(245) nov.(214) year(208) at(183) sept.(179) #(163) stake(158) outstanding(151) common(150) mark(139) ton(129) one(125) <unk>(119) franc(119) nearly(112) march(112) least(110) dec.(107) off(105) fiscal(103) each(101) only(86) three(86) dollar(85) five(81) c(80) pence(79) record(78) barrel(76) july(74) in(71) April(69) six(68) per(66) chapter(64) between(62) p.m.(62) a.m.(60) due(59) percentage(59) trillion(58) no.(57) aug.(57) additional(56) series(55) foot(55) june(54) Canadian(54) two(51) mile(51) accord(49) around(48) fall(47) jan.(47) October(46) ounce(46) June(45) bid(44) late(44) four(43) rise(43) annually(43) pound(42) almost(41) last(40) an(39) roughly(39) Swiss(39) just(37) hour(37) metric(37) minute(35) age(35) may(34) add(34) page(33) week(32) people(32) Tuesday(32) fix(31) basis(30) nyse(29) be(28) unit(28) equal(28) warrant(28) which(27) 
[7]
<unk>(9382) new(2510) ##(1277) u.s.(941) first(825) stock(755) big(716) federal(621) major(551) two(547) other(530) next(511) financial(416) national(380) good(372) last(364) past(361) most(359) American(358) recent(353) few(345) large(345) own(335) third(325) Japanese(307) real(306) investment(305) small(303) economic(301) same(298) san(293) late(287) business(284) current(274) market(266) second(259) high(257) trade(250) soviet(242) tax(241) security(240) state(233) British(231) wall(224) more(218) three(218) future(217) foreign(207) nine(201) public(199) computer(198) political(196) white(192) international(192) junk(190) trading(178) previous(174) insurance(172) oil(172) early(169) great(165) European(165) takeover(163) old(162) exchange(160) corporate(159) fiscal(154) world(152) california(152) strong(151) government(150) money(149) only(149) auto(148) los(146) very(146) german(145) former(141) low(137) home(137) legal(136) private(136) propose(135) house(133) bond(132) joint(131) holding(131) special(130) five(129) capital(129) program(129) west(126) long(125) bank(125) defense(124) local(123) treasury(123) common(122) credit(120) work(119) fourth(118) 
[8]
'(12315) s(6104) S(4711) re(389) ve(200) jones(194) m(143) ll(130) industrial(122) d(104) &(60) say(56) lambert(49) hutton(48) burnham(43) lehman(36) /(27) Lehman(23) capital(19) ua(15) airline(13) dow(13) lynch(12) <unk>(11) equity(11) professional(11) Co(10) Hutton(10) estimate(9) .(8) ready(8) earlier(8) expectation(8) u(8) transportation(7) container(7) Burnham(7) 2-year(7) manager(6) inc(6) analyst(6) brother(6) bradstreet(6) firm(5) line(5) financial(5) fund(5) investment(5) saving(5) dealer(5) home(5) marketing(4) report(4) share(4) union(4) decline(4) security(4) add(4) lawyer(4) committee(4) index(4) retain(4) organization(4) compensation(4) casualty(4) new(3) right(3) mr.(3) producer(3) group(3) trading(3) meeting(3) machinist(3) average(3) agency(3) air(3) product(3) future(3) loan(3) e(3) regulator(3) action(3) player(3) worker(3) manufacturer(3) science(3) social(3) instrument(3) dataproducts(3) brand(3) calif(3) symptom(3) cruise(3) resident(3) to(2) for(2) ago(2) start(2) park(2) phone(2) fee(2) 
[9]
<unk>(6878) ##(1482) more(1100) it(1053) out(870) up(789) them(689) much(579) one(564) all(564) such(466) well(464) some(398) accord(397) back(392) sale(382) because(374) him(349) part(330) most(283) down(272) many(270) stock(261) price(236) which(229) order(215) yesterday(209) money(201) time(199) addition(198) control(197) and(196) what(196) only(195) example(192) congress(190) September(190) japan(189) people(188) investor(188) that(186) those(181) august(179) how(179) us(176) as(175) now(171) less(171) business(170) off(169) least(163) even(159) today(156) cash(154) just(151) away(151) interest(147) good(146) Friday(143) bank(141) base(139) profit(137) trading(137) asset(136) half(136) work(134) this(129) high(128) europe(126) again(125) not(123) concern(120) comment(120) other(118) charge(118) year(115) trade(115) here(115) itself(115) home(114) talk(111) debt(111) demand(110) Monday(109) soon(109) ahead(109) me(109) term(108) early(107) long(106) share(106) two(105) so(105) loan(105) on(104) damage(102) earning(102) rather(102) place(98) along(97) change(96) 
[10]
<unk>(5803) year(1296) market(700) company(583) time(514) month(496) week(480) sale(441) day(415) price(377) number(369) way(359) quarter(351) share(341) stock(336) result(332) cost(329) end(326) trading(325) u.s.(323) loss(288) case(279) issue(275) problem(266) plan(262) business(255) group(255) board(243) bank(235) office(234) level(232) right(226) offer(226) rate(225) value(220) bid(217) lot(216) agreement(216) gain(216) contract(214) position(207) unit(200) effort(196) security(193) interest(183) product(183) state(180) investment(178) increase(176) house(169) decline(168) operation(168) government(165) bill(163) amount(161) maker(158) stake(157) firm(155) meeting(154) decision(152) system(151) change(151) part(151) program(150) point(150) country(148) use(145) proposal(144) economy(144) area(143) report(141) member(137) line(136) plant(134) period(133) fund(133) deal(131) policy(130) move(129) thing(128) director(126) kind(125) charge(124) president(123) series(123) acquisition(122) good(119) purchase(117) department(117) hour(117) asset(116) risk(116) return(114) court(114) loan(114) role(114) average(111) effect(111) debt(110) job(110) profit(109) 
[11]
be(17264) have(6922) will(3396) n't(3342) would(2471) which(1863) do(1809) who(1650) that(1520) also(1504) could(1210) can(1016) may(888) not(736) and(571) should(471) still(452) now(402) might(383) <unk>(375) already(307) say(284) must(269) wo(263) ca(206) official(198) just(188) help(186) never(186) currently(158) recently(157) even(147) often(146) probably(145) corp.(138) then(137) only(120) he(109) to(108) begin(106) really(105) remain(104) Inc.(96) yet(95) actually(92) always(88) generally(85) it(84) however(84) all(81) ago(80) simply(78) they(75) ever(75) rate(74) long(73) apparently(71) usually(69) but(68) previously(68) old(68) become(65) well(65) once(62) seem(61) here(57) eventually(56) certainly(56) soon(53) typically(53) further(51) later(51) Corp.(51) spokesman(50) group(49) start(47) quickly(47) consider(46) finally(46) first(45) sometimes(45) itself(44) before(43) order(43) you(43) either(43) report(42) feel(42) i(41) so(40) executive(39) holder(39) too(39) without(39) fully(38) immediately(38) clearly(37) after(36) stop(35) show(34) more(33) 
[12]
it(5422) he(3796) <unk>(2759) they(2657) we(1494) i(1289) there(1022) that(780) you(735) she(555) this(453) analyst(369) people(271) bush(236) investor(233) trader(226) what(203) japan(181) yesterday(179) Friday(143) congress(115) ibm(112) moody(109) ford(93) china(92) Drexel(90) warner(90) today(88) dealer(80) other(77) peter(76) Noriega(75) dow(74) britain(72) gm(72) not(71) ual(69) those(67) dinkins(66) one(64) Gorbachev(62) official(61) here(61) some(60) jaguar(59) price(57) Monday(56) who(55) thatcher(54) Tuesday(53) term(50) columbia(50) sony(49) Lloyd(49) American(47) why(47) Nissan(47) france(47) economist(47) ##(46) London(46) baker(45) eastern(45) Shearson(43) these(42) Paul(42) digital(42) cbs(40) after(40) shearson(40) everyone(39) Phelan(39) california(38) sear(38) australia(37) most(36) Sotheby(36) trump(36) kidder(36) many(35) nobody(34) Krenz(34) consumer(33) bank(33) someone(33) both(32) rally(32) usx(32) smith(32) Lawson(32) but(31) anyone(30) Canada(30) mccaw(30) lawmaker(30) lang(30) do(29) trading(29) jones(29) then(29) Exxon(29) 
[13]
<unk>(6684) be(5433) have(2372) make(1590) sell(1088) take(1050) expect(956) get(929) go(846) buy(766) use(713) pay(660) come(570) do(567) want(555) continue(549) give(533) not(510) see(504) hold(438) plan(423) call(416) try(390) begin(386) provide(383) offer(378) put(371) agree(368) raise(357) become(347) work(344) receive(338) seek(336) increase(332) need(327) n't(319) find(313) keep(312) report(310) run(301) no(299) own(299) acquire(297) help(293) in(291) look(286) include(286) build(282) reduce(277) close(275) turn(275) set(274) show(272) move(270) consider(269) decline(267) lead(264) allow(263) reach(260) know(258) say(247) produce(247) leave(246) likely(245) file(243) so(237) remain(236) cut(236) lose(235) create(233) meet(232) tell(231) seem(230) start(226) require(226) bring(225) spend(220) more(220) open(211) fall(210) add(207) able(206) appear(201) end(196) announce(195) force(194) win(191) change(191) develop(190) only(189) like(189) cause(187) complete(187) just(181) approve(178) very(175) grow(172) ask(172) face(172) follow(167) base(165) 
[14]
the(54731) a(20191) its(4115) an(3347) their(1960) his(1948) this(1944) <unk>(1941) some(1153) other(1054) last(811) any(801) ##(760) that(754) one(734) no(589) these(574) many(542) two(519) our(445) another(412) such(394) her(386) new(384) those(380) several(357) more(351) all(327) three(319) u.s.(318) each(273) both(265) my(265) most(258) federal(236) your(235) mr.(226) every(215) high(214) foreign(188) certain(177) program(175) four(160) American(160) five(154) recent(154) big(147) Japanese(144) early(140) small(125) major(124) six(116) low(116) next(115) large(115) south(110) wall(103) market(101) commercial(98) corporate(98) future(97) west(95) national(94) stock-index(93) east(90) British(89) additional(81) strong(80) which(75) interest(75) state(75) international(74) western(74) great(72) good(72) further(72) long-term(72) short-term(70) composite(69) little(68) various(68) seven(67) local(66) government(66) very(64) increase(63) late(63) consumer(61) hong(61) least(60) economic(60) individual(59) it(58) potential(58) general(57) improve(56) private(56) full(53) index(52) personal(51) eight(49) 
[15]
<unk>(3617) company(2968) year(2156) market(1525) york(1004) bank(631) group(593) firm(588) month(577) business(519) fund(512) government(511) quarter(508) industry(490) week(475) investor(443) price(441) stock(431) board(424) issue(414) product(402) and(396) exchange(369) official(342) agency(342) bill(336) plan(317) system(308) court(308) bond(303) time(302) program(299) house(297) state(281) trading(267) concern(265) street(265) department(264) unit(263) country(263) nation(259) service(256) operation(247) u.s.(241) people(237) administration(235) union(228) day(227) sale(226) maker(217) problem(215) contract(212) rate(211) committee(210) share(209) law(206) plant(202) dollar(198) commission(198) period(192) spokesman(191) francisco(191) estate(191) analyst(187) area(186) economy(182) report(173) city(173) world(172) trader(172) dow(172) agreement(167) one(163) manager(163) airline(159) case(154) move(153) loan(153) man(153) project(149) policy(149) offer(147) computer(144) transaction(144) worker(143) judge(139) line(137) bid(136) debt(134) executive(133) party(131) security(131) deal(131) p(129) change(128) figure(125) woman(125) treasury(124) angeles(123) leader(123) proposal(123) 
[16]
to(16906) and(2332) or(1621) from(1347) by(229) month(183) but(163) us(88) longer(60) than(40) high(37) without(37) note(31) up(30) while(30) about(26) low(26) so(26) off(20) turn(18) britain(15) what(14) utility(14) japan(13) germany(12) Switzerland(12) transportation(12) difficulty(12) far(11) out(10) problem(10) fully(10) order(9) set(8) phone(7) before(7) you(7) against(7) sheet(7) payable(7) them(6) less(6) effectively(6) wsj(6) its(5) another(5) plunge(5) around(5) pay(5) back(5) instruction(5) me(5) bad(5) him(5) produce(5) their(4) this(4) since(4) amid(4) publicly(4) thus(4) product(4) fairly(4) closely(4) overseas(4) thereby(4) telephone(3) start(3) people(3) more(3) see(3) further(3) activity(3) virtually(3) upon(3) magazine(3) factor(3) economic(3) themselves(3) sometimes(3) item(3) formerly(3) fashion(3) room(3) slowly(3) traditionally(3) pollution(3) whichever(3) consumer(2) can(2) record(2) technology(2) with(2) way(2) sell(2) ad(2) system(2) director(2) if(2) customer(2) toward(2) 
```

```
alpha 9.9941e-05
beta[1] 1.81592
beta[2] 0.216374
beta[3] 2.21947
beta[4] 5.14212
beta[5] 2.68645
beta[6] 0.983193
beta[7] 3.23531
beta[8] 0.961603
beta[9] 4.40115
beta[10] 3.11974
beta[11] 4.09954
beta[12] 1.87855
beta[13] 1.33697
beta[14] 3.63117
beta[15] 5.71038
beta[16] 0.318908
```

### こゝろ

`kokoro.txt`の学習結果です。

```
python3 train_ja.py  -f ../../text/kokoro.txt -split 1 --start-temperature 1.5 --min-temperature 0.08 -tags 16 -alpha 1 -beta 1 -e 40000
```

![result](https://raw.githubusercontent.com/musyoku/images/master/pos_tagging/bhmm/kokoro.png)

```
[1]
て(3323) ば(200) た(109) ながら(83) たら(67) たり(59) なく(52) ない(50) と(26) たく(21) で(20) ので(20) ず(14) ませ(13) られ(12) ちゃ(11) <unk>(9) かも(8) つつ(8) って(7) し(6) ら(5) のに(4) だら(4) らしく(3) 得(3) か(3) 得る(3) れる(3) たい(3) おっしゃる(3) も(2) さ(2) け(2) なる(2) 間(2) とか(2) と共に(2) べき(2) てる(2) られる(2) 気遣い(2) 更(2) へ(1) う(1) かえって(1) あり(1) 自由(1) 以上(1) なけれ(1) そう(1) 
[2]
<unk>(891) 私(691) 先生(261) それ(244) 自分(225) Ｋ(202) 何(177) 彼(173) 人(150) 奥さん(141) 父(105) そこ(103) お嬢さん(91) あなた(85) 母(83) ない(73) 東京(69) そう(62) 口(60) 今(57) 心(53) 急(52) 気(51) 宅(47) 眼(43) 手(42) 妻(41) どこ(38) 叔父(37) 外(35) 変(35) 女(35) いっしょ(35) 前(33) 他(32) ここ(31) 頭(31) いつも(30) 世の中(30) 少し(30) 話(29) 先(29) 的(28) 時(28) 仕方(28) 向う(28) ##(27) 」(27) 顔(26) 日(26) 手紙(26) 
[3]
私(1915) 先生(299) 奥さん(224) Ｋ(194) 父(139) それ(135) 彼(129) あなた(79) 時(69) <unk>(66) 人(66) 母(63) お嬢さん(62) これ(60) 妻(53) 兄(34) 今(28) 我々(27) 自分(21) 君(20) お前(17) 叔父(17) 医者(17) 頃(17) 今度(16) 日(16) おれ(16) 家(13) 人間(13) 女(13) そこ(12) 手紙(12) 眼(11) 言葉(11) もの(11) 彼ら(11) 「(10) 何(10) どこ(10) 時分(10) 誰(10) お父さん(10) ##(8) 話(8) 恋(8) 男(8) しまいに(8) 子供(8) 中(7) 他(7) 態度(7) 
[4]
た(5213) です(1087) ん(383) ない(357) う(266) ます(238) だ(195) ね(66) よ(58) たい(26) られる(23) ある(22) か(21) 下さい(19) から(15) ちゃ(14) な(13) …(12) さ(11) わ(11) いる(9) まい(9) なる(8) れる(8) てる(8) する(7) ば(7) たろ(5) なさい(5) <unk>(4) て(4) ）(4) かい(4) 廻る(4) は(3) せる(3) ず(3) すか(3) 出す(3) みえる(3) も(2) 」(2) 過ぎる(2) って(2) 上げ(2) べき(2) 付く(2) 頂戴(2) ください(2) たまえ(2) させよ(2) 
[5]
は(3768) も(607) が(452) の(138) から(87) と(72) 、(71) ##(66) を(34) か(29) に(27) なく(26) でも(23) なら(20) に対する(17) だって(14) や(13) ない(12) まで(12) う(9) に対して(9) とも(9) <unk>(8) ので(8) ば(8) まま(8) って(7) …(7) で(6) また(6) 以来(6) 人(4) 自身(4) ね(4) 時(4) どうして(4) に関して(4) かも(3) です(3) じゃ(3) よく(3) どう(3) 少し(3) ごとく(3) 大変(3) もっと(3) 手前(3) のみ(3) た(2) し(2) あり(2) 
[6]
<unk>(874) いる(336) ない(232) ある(160) する(138) いう(118) も(99) なる(82) 同じ(71) ##(68) 来る(59) どう(50) 思う(47) 行く(42) よく(40) 見る(40) いわ(39) 好い(37) いい(37) 卒業(33) 思わ(32) その(32) どこ(30) また(28) そう(28) 出る(28) 悪い(27) そんな(27) さ(25) 聞か(25) 聞く(24) 帰る(24) こう(24) この(23) 強い(22) まるで(21) 何(21) みる(21) 死ぬ(20) と(19) 私(19) 想像(18) し(17) 落ち(17) まだ(17) こんな(17) くれる(17) 人(16) 通り(16) それ(15) という(15) 
[7]
、(3416) も(114) いる(83) ##(78) また(55) その(51) <unk>(27) よく(22) まだ(20) この(19) から(17) ある(17) ほとんど(15) 少し(12) ただ(11) すぐ(10) お(10) たった(10) ない(9) かえって(9) なら(9) もう(9) 小(9) する(8) 「(7) さえ(7) むしろ(7) くれる(7) 新しい(6) 決して(6) わざと(6) 常に(6) なく(5) 無論(5) 突然(5) 大きな(5) わざわざ(5) 行く(5) 実際(5) しか(5) きっと(5) よほど(5) とうとう(5) こう(5) それほど(5) と(4) は(4) まるで(4) みんな(4) 今(4) 一つ(4) 
[8]
し(1132) <unk>(538) あっ(321) いっ(259) れ(234) なっ(213) 見(128) 出(119) 思っ(108) でき(100) 聞い(100) 考え(96) なかっ(93) なけれ(88) 帰っ(87) い(81) 来(80) 出し(69) 見え(67) あり(66) 答え(60) 行っ(56) もっ(55) られ(53) せ(48) 立っ(48) 書い(45) 話し(41) 坐っ(39) ない(37) 知っ(37) 違っ(34) 持っ(34) も(32) 黙っ(32) 知れ(31) いい(31) 向っ(31) 信じ(31) 付い(30) 眺め(29) 感じ(29) いえ(27) 取っ(27) 生き(26) 思い(26) 開け(26) 始め(25) ん(24) 置い(24) 尋ね(24) 
[9]
い(766) あり(157) 来(154) くれ(72) なり(57) いい(56) なら(56) 行っ(46) ん(45) いる(44) しまっ(44) なかっ(43) 知れ(40) 思い(35) なっ(33) <unk>(32) み(32) あっ(28) しまい(24) し(20) 見(19) 聞き(18) 出(17) おい(17) いけ(16) 見せ(15) おき(13) やり(13) 行き(12) もらっ(12) いっ(10) 思わ(10) 解ら(10) もの(10) 答え(10) 見え(10) 笑い(9) なれ(9) 聞い(9) 驚き(9) 考え(8) 思え(8) き(8) 同じ(8) 分り(8) なかろ(8) 始め(7) 違い(7) 起ら(7) も(6) 帰っ(6) 
[10]
まし(1030) ませ(353) なかっ(327) でし(222) だっ(43) でしょ(42) られ(40) ん(39) 出し(36) て(29) れ(28) ましょ(28) だろ(26) らしかっ(23) たかっ(19) だ(16) です(13) 始め(12) なくっ(11) そう(11) 事(8) せ(8) ます(6) たく(6) <unk>(4) 結構(4) 過ぎ(4) 聞こえ(4) っ(4) 浮べ(4) ない(3) 。(3) 得(3) よかっ(3) ござん(3) の(2) し(2) たろ(2) は(2) よう(2) 帰っ(2) か(2) なり(2) じゃ(2) つもり(2) 直し(2) つけ(2) 込ん(2) 繰り返し(2) 頃(2) てる(2) 
[11]
の(5186) な(585) という(310) ない(134) でし(114) だろ(113) でしょ(105) する(99) から(79) ん(77) だ(65) で(46) もの(46) た(41) より(40) が(37) ##(36) れる(36) らしい(35) べき(35) か(31) といった(31) と(30) や(28) くらい(21) せる(19) る(18) よう(17) 事(16) に(15) まで(15) じゃ(15) に対する(15) がる(15) <unk>(14) その(11) い(9) でも(9) はず(9) とかいう(9) ぎり(9) とも(9) たい(8) しよ(8) です(7) も(7) 」(7) ほど(7) られる(7) れ(6) 大きな(6) 
[12]
に(4035) を(3172) と(1816) が(1692) で(1418) から(580) へ(512) も(432) は(358) の(253) か(167) まで(122) より(108) として(96) でも(94) な(86) 、(85) じゃ(65) について(61) さ(59) <unk>(57) かも(51) し(45) なら(43) だけ(42) さえ(41) や(40) に対して(36) ほど(34) 通り(28) だ(25) ばかり(25) のに(20) とも(20) だの(16) らしく(15) 時(15) にとって(15) とか(14) たち(13) しか(13) って(13) なく(12) たり(12) と共に(11) において(11) ので(10) 」(9) ば(9) なり(9) ず(9) 
[13]
「(677) その(586) ##(315) しかし(210) この(182) また(167) <unk>(146) そうして(125) ただ(98) けれども(77) 時(76) すぐ(73) まだ(61) すると(60) もう(56) それから(52) それで(41) 、(40) 今(40) だから(39) そう(39) そんな(33) もし(32) むしろ(32) あの(31) 突然(31) 日(31) それでも(30) こう(30) ある(29) こんな(29) 実際(29) ――(27) なぜ(26) 上(26) かえって(23) もっとも(23) 無論(23) 全く(23) 時々(23) ことに(22) 決して(22) あるいは(22) ところが(22) ちょっと(21) これから(20) こういう(20) お(20) 若い(19) とうとう(19) つまり(19) 
[14]
<unk>(271) 事(167) よう(150) 方(142) もの(122) 中(107) 前(99) 顔(97) 上(80) うち(78) ため(76) 言葉(71) 間(60) 人(57) 眼(53) う(52) 室(50) 心(49) 私(39) 態度(38) 意味(34) 頭(34) 様子(33) 家(32) 声(32) 所(31) の(29) 口(28) 胸(28) 病気(28) 過去(26) 」(25) 心持(24) 気(21) 手(21) だ(20) 手紙(20) 話(20) 父(19) 外(19) 返事(19) 調子(19) 気分(19) 色(19) 通り(18) ところ(17) 希望(17) 後(16) 下(16) 名(16) 傍(16) 
[15]
。(4632) 」(521) か(186) 時(83) から(57) ので(53) けれども(52) よ(41) ね(28) 、(15) まま(12) 後(10) のに(9) もの(8) が(7) らしい(7) なら(5) ）(4) かい(4) わ(4) ため(3) その(2) 通り(2) ぜ(2) …(1) がっ(1) 堅く(1) 市(1) 晩(1) 本当に(1) だい(1) とき(1) 末(1) 
[16]
事(349) よう(326) だ(262) か(202) う(201) もの(179) の(178) <unk>(171) だけ(103) 人(87) ない(86) 時(83) ず(73) 気(64) 」(54) 方(43) ところ(43) と(41) 訳(34) など(33) ため(32) 通り(31) 所(30) さ(30) くらい(30) ほど(30) ばかり(29) する(28) 後(28) 男(27) うち(26) つもり(26) 前(25) そう(22) はず(21) 家(19) 私(19) 人間(19) れる(19) 的(18) まで(18) 日(18) 機会(18) 様子(17) 意味(17) 風(17) た(16) 必要(16) 問題(16) 心持(16) 言葉(15)
```

```
alpha 9.45495e-08
beta[1] 0.0201575
beta[2] 10.5256
beta[3] 3.56328
beta[4] 0.0216688
beta[5] 0.0130361
beta[6] 7.25476
beta[7] 2.20238
beta[8] 4.70198
beta[9] 0.0964968
beta[10] 0.161765
beta[11] 0.694262
beta[12] 0.0140772
beta[13] 5.41368
beta[14] 4.76062
beta[15] 0.0665414
beta[16] 4.30079
```

### 吾輩は猫である

`neko.txt`の学習結果です。

```
python3 train_ja.py  -f ../../text/neko.txt -split 1 --start-temperature 1.5 --min-temperature 0.08 -tags 16
```

![result](https://raw.githubusercontent.com/musyoku/images/master/pos_tagging/bhmm/neko.png)

```
[1]
主人(619) <unk>(510) 吾輩(316) これ(248) 君(248) それ(236) 何(196) 僕(172) 迷亭(144) 細君(144) 人間(132) 人(95) 自分(91) 私(89) 今(81) 彼(80) 例(67) 今度(65) 独(62) 苦(61) 猫(60) 誰(60) 先生(57) 昔(57) 今日(51) あれ(48) あなた(40) 彼等(34) 寒月(34) ここ(33) どこ(33) そこ(32) 右(31) 鼻(31) 妙(31) かく(29) 金田(29) ところ(28) うち(28) 云う(27) 仕方(27) 女(26) 敵(26) 見る(24) 気(23) こっち(23) 学校(22) 御前(21) 少し(21) 世の中(21) ##(20) 
[2]
<unk>(499) い(274) 来(174) も(164) 見(136) しまっ(78) なら(65) くれ(49) し(47) は(43) どう(40) 行っ(38) 出来(32) やり(31) 思わ(30) と(27) 云わ(27) 出(26) ちょっと(26) 行か(26) おい(25) よく(25) なく(22) おっ(21) 極(21) おら(21) やっ(20) でも(20) なり(17) どうか(16) まだ(16) 相(16) 発達(16) にやにや(16) なっ(15) おり(15) 感心(15) 聞か(15) 知ら(14) そう(14) 早く(14) やら(14) あっ(13) 聞い(13) こう(12) しまい(12) 一つ(12) 拝見(12) しばらく(12) のみ(12) 取ら(12) 
[3]
<unk>(613) 君(400) 事(369) 人(270) さん(215) ##(161) もの(130) 時(129) よう(120) だ(116) など(109) の(103) 子(88) くらい(87) だけ(85) そう(76) 声(74) 者(72) 家(69) 気(67) 日(66) 方(65) 先生(64) ところ(51) 度(51) 年(49) 様(49) 供(49) 顔(47) か(45) 」(42) 事件(39) 所(39) 男(39) 本(39) まで(38) 内(38) 的(37) 屋(37) 館(34) 中(33) 猫(32) 通り(32) 金(31) 訳(31) 頭(31) ず(30) 目(30) 奴(29) 女(29) 以上(28) 
[4]
し(1902) <unk>(908) なっ(363) れ(285) 云っ(225) 思っ(211) 出(192) あっ(189) まし(189) 見(176) 来(162) 聞い(125) 見え(117) 出し(117) 出来(115) い(112) せ(105) なけれ(102) 行っ(98) なく(98) なかっ(96) 持っ(92) られ(79) もっ(77) なら(74) 考え(69) つけ(63) やっ(59) すれ(59) 這入っ(57) 知ら(56) 食っ(56) 入れ(55) っ(54) 得(51) 知っ(51) 云わ(50) なくっ(47) 寝(44) 帰っ(44) 相違(42) 立っ(41) 云い(37) 上っ(37) 始め(37) 買っ(35) 笑っ(35) かけ(35) かい(34) 云え(33) 致し(30) 
[5]
<unk>(636) 云う(273) する(247) なる(201) ある(150) 思う(66) 見る(61) 聞く(57) 見える(49) 出す(42) 出る(34) 主人(30) 何(30) ない(28) いう(25) 行く(24) やる(24) 人(22) 吾輩(22) ええ(21) いい(20) 君(19) 心(19) 出来る(19) 馬鹿(18) 相手(17) ##(17) よる(17) 読ん(16) 僕(15) 廻る(15) 猫(14) す(14) 手(14) くる(14) 細君(14) なるほど(14) これ(13) かける(13) 飲ん(13) 気(13) 帰る(13) ここ(12) それ(12) 食う(12) 書斎(11) 顔(11) 答える(11) 入る(11) 尻(11) 引く(11) 
[6]
<unk>(1138) 事(262) よう(230) 上(186) 方(163) 中(149) 顔(125) 何(122) もの(104) 主人(104) 人(92) 頭(84) ない(81) 君(76) 前(73) ため(70) うち(69) ところ(68) 人間(64) 大(64) 気(64) これ(61) 猫(57) 眼(57) の(54) 手(54) 声(52) 間(50) 口(48) 所(47) 鼻(43) 急(41) どこ(40) 下(40) 吾輩(39) 家(37) ##(34) 先生(32) 音(32) 」(31) ある(31) ほか(31) 先(31) 裏(29) ここ(27) 自分(27) 横(27) それ(26) 細君(26) ヴァイオリン(26) 名(26) 
[7]
だ(1335) ある(1010) ない(830) う(704) か(604) ね(578) いる(574) よ(459) です(449) ん(330) さ(304) な(190) …(155) かい(99) た(90) ます(89) ぜ(84) ？(79) いい(77) わ(67) だい(57) から(49) くる(40) ねえ(39) もの(37) の(32) なる(30) さあ(30) 来る(29) なさい(27) まい(26) <unk>(25) 見る(23) しまう(20) 下さい(20) は(18) 行く(18) って(18) らしい(17) え(16) 面白い(15) あ(15) なあ(15) で(14) ぬ(14) 云う(13) や(12) 君(12) する(11) なり(11) 困る(11) 
[8]
を(5324) に(4996) と(2119) が(1716) で(1473) へ(885) から(497) は(423) も(414) まで(141) の(136) か(108) さ(87) として(82) さえ(73) ばかり(69) だけ(65) し(63) でも(61) て(43) ほど(42) ろ(42) なく(39) とか(39) や(32) なんか(32) <unk>(31) より(30) 、(28) など(22) くらい(21) なら(19) 々(19) 暗に(18) によって(17) り(16) だって(16) 人(15) ども(15) よく(15) ら(15) って(14) ず(13) ながら(13) において(13) れ(12) やら(12) せら(12) ば(11) 度(11) のみ(11) 
[9]
の(7717) は(4178) が(2342) な(1227) も(896) で(665) に(601) から(519) と(332) か(280) 、(279) なら(202) でも(191) だ(184) を(137) じゃ(132) より(124) なる(118) や(108) する(98) ##(92) <unk>(90) ば(88) だって(84) まで(81) 沙弥(76) という(73) 仙(72) ので(64) たる(56) ばかり(52) し(50) って(49) ながら(41) 的(37) として(36) だけ(33) において(32) べき(31) ほど(31) とか(31) 衛門(31) へ(28) のに(28) に対して(24) にゃ(23) について(22) くらい(20) 共(20) ぬ(17) る(17) 
[10]
「(2711) しかし(114) ――(93) この(63) ただ(54) すると(45) その(42) ##(37) だから(34) こう(32) もし(31) 第(24) それから(23) ところが(23) 元来(21) いくら(20) また(19) それで(18) もっとも(18) なるほど(17) 寒月(17) 今(16) そこで(16) 実は(15) もう(14) どうも(14) そうして(14) 現に(14) ことに(13) しかも(13) ねえ(12) まあ(11) 但し(11) 君(10) 御(10) あの(10) …(10) しかるに(10) 何だか(9) いかに(9) 東風(9) せんだって(9) どうせ(9) なぜ(9) 従って(9) と(8) 先生(8) ああ(8) どう(8) やはり(8) ちょうど(8) 
[11]
で(1176) ん(799) だ(481) です(416) もの(304) が(297) ませ(278) は(265) じゃ(262) ない(253) だろ(245) そう(161) も(161) …(154) でしょ(152) 事(135) 知れ(127) かも(107) あり(105) から(96) か(95) よう(92) <unk>(77) の(76) あろ(67) もん(62) ある(59) う(58) ござい(54) なら(48) 出来(45) いい(45) どう(45) ます(43) 男(42) 分ら(41) 訳(41) ましょ(40) はず(38) ところ(34) 駄目(34) よかろ(32) ばかり(28) いけ(28) くらい(26) つもり(25) いか(24) 知ら(23) でも(21) あ(21) なかろ(20) 
[12]
##(1035) <unk>(906) 云う(475) この(435) 御(421) その(408) そう(179) そんな(172) 寒月(166) あの(153) ない(128) ある(127) する(101) こんな(96) なる(93) 君(93) 迷亭(90) また(85) 大(74) 大きな(73) 東風(72) まあ(72) 小(72) ええ(72) まだ(71) 鈴木(67) ちょっと(66) 少し(66) 金田(66) もう(62) いい(60) どう(59) そりゃ(54) 先生(50) どうも(50) 時(48) よく(47) 第(46) 鼻(45) しかし(45) 同じ(44) いや(43) なるほど(43) それから(43) そんなに(42) ただ(42) それで(42) あんな(42) 何(40) なに(40) なかなか(40) 
[13]
。(6811) 」(2895) から(558) が(420) か(133) ので(46) まい(31) って(27) のに(27) けれども(26) と(20) とも(16) など(15) もの(11) で(10) …(8) らしい(8) 、(7) なら(6) ね(6) に(5) ものの(5) けれど(5) ――(4) よ(4) ぜ(4) なんて(4) <unk>(3) 間(3) て(3) よう(3) い(3) 訳(3) 時分(3) んで(3) つもり(3) は(2) の(2) うち(2) し(2) 退校(2) え(2) わす(2) 景色(2) 騒々しい(2) ながら(1) 事(1) だ(1) さ(1) れる(1) ん(1) 
[14]
いる(375) もの(304) 事(293) <unk>(289) の(186) よう(176) ある(127) ない(125) 見る(103) あり(75) ところ(70) 主人(66) いい(66) 時(66) か(59) くらい(38) 者(38) 訳(38) 吾輩(34) 方(33) 男(32) 人(30) 何(30) ばかり(26) 仕方(25) くる(25) 顔(23) 出来る(22) ん(20) うち(19) 日(19) あと(19) 猫(18) やる(18) くれる(18) 間(17) おく(17) くれ(16) 所(16) つもり(16) しまう(15) 人間(15) 上(15) 以上(15) 眼(14) これ(13) 行く(13) 鼻(13) 通り(13) ため(13) 例(13) 
[15]
て(6202) た(3574) ない(716) ば(448) する(359) ます(258) ん(238) たら(233) ず(205) ぬ(205) てる(160) う(128) たい(109) ながら(102) たり(102) ちゃ(80) れる(73) ざる(66) まし(44) 候(39) せる(37) と(35) る(35) だ(34) たく(34) たる(33) るる(31) <unk>(30) たって(30) の(29) ませ(27) なく(25) ましょ(23) ら(23) に(22) す(22) なる(21) 給え(21) られ(20) で(19) ねえ(18) そう(17) って(17) 得る(17) つ(17) つつ(17) たまえ(17) たろ(16) なさい(15) られる(14) な(12) 
[16]
、(5887) と(2373) は(809) も(598) ##(394) 「(263) <unk>(207) ――(177) いる(164) ごとく(148) から(127) 御(116) また(95) この(70) って(63) が(62) その(50) か(45) ちょっと(43) まま(42) 時(41) 必ず(39) ある(37) ので(36) なく(35) しきりに(28) に(27) もう(27) 少々(27) ？(26) 第(26) よく(26) やはり(25) 少し(24) 決して(23) 通り(23) お(23) 迷亭(22) する(20) 云う(19) あまり(18) 寒月(18) 」(17) を(17) で(17) 到底(17) くらい(17) 見る(16) 何だか(15) 突然(15) 皆(15) 
```

```
alpha 4.57117e-09
beta[1] 6.79579
beta[2] 4.96007
beta[3] 6.28583
beta[4] 3.67425
beta[5] 6.46577
beta[6] 8.72178
beta[7] 0.100578
beta[8] 0.0356582
beta[9] 0.221254
beta[10] 0.387118
beta[11] 1.47342
beta[12] 7.06867
beta[13] 0.0261711
beta[14] 3.86126
beta[15] 0.0614419
beta[16] 3.60394
```