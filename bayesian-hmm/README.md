## A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging

- [A Fully Bayesian Approach to Unsupervised Part-of-Speech Tagging](http://homepages.inf.ed.ac.uk/sgwater/papers/acl07-bhmm.pdf)
- [実装について](http://musyoku.github.io/2017/01/28/A-Fully-Bayesian-Approach-to-Unsupervised-Part-of-Speech-Tagging/)

## 準備

### macOS

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
./bootstrap.sh --with-libraries=python --with-python=python --with-python-root=YOUR_PYTHON_ROOT
./b2 python=2.7 -d2 -j4
./bootstrap.sh --with-libraries=python --with-python=python3 --with-python-root=YOUR_PYTHON_ROOT
./b2 python=3.6 -d2 -j4
```

```
make install
```

### MeCabのインストール

```
pip install mecab-python3
```

## 学習

### 英語

```
python train_en.py -f ../alice.txt -n 7
```

### 日本語

```
python train_ja.py -f ../aozora.txt -n 20
```

### 結果の可視化

```
python tags.py -n 200
```

獲得した品詞とそれに属する単語を一覧表示します。

## 結果

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