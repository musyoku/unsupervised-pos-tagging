## The Infinite Hidden Markov Model

- [The Infinite Hidden Markov Model](http://mlg.eng.cam.ac.uk/zoubin/papers/ihmm.pdf)
- [Beam Sampling for the Infinite Hidden Markov Model](http://mlg.eng.cam.ac.uk/pub/pdf/VanSaaTehGha08.pdf)
- [実装について](http://musyoku.github.io/2017/02/27/Infinite-Hidden-Markov-Model%E3%81%AB%E3%82%88%E3%82%8B%E6%95%99%E5%B8%AB%E3%81%AA%E3%81%97%E5%93%81%E8%A9%9E%E6%8E%A8%E5%AE%9A/)

#### Todo:

- [ ] 日本語学習用のPythonコード
- [ ] ハイパーパラメータのサンプリング
- [ ] Viterbiデコード
- [ ] ビームサンプリング

## ビルド

```
make install
```

## 学習

### 英語

```
python train_en.py -f ../alice.txt -n 2 -e 100000 -u 1
```

## 結果

```
tag 1:
	<eos>/1397, 
tag 2:
	./1209, !/115, ?/74, 
tag 3:
	queen/75, king/64, head/60, hatter/57, gryphon/55, mouse/48, duchess/42, dormouse/40, eye/36, door/32, moment/31, caterpillar/27, hand/26, jury/22, arm/21, other/21, face/20, word/20, house/19, cat/19, court/18, table/18, air/16, tree/15, footman/15, baby/14, sea/14, game/14, dance/14, mouth/14, lobster/13, life/13, cook/13, dodo/13, tail/12, dear/12, majesty/12, slate/12, pigeon/12, name/12, trial/11, bottle/10, hedgehog/10, soldier/10, reason/10, pool/10, rate/9, knave/9, history/8, direction/8, window/8, shoulder/8, distance/8, wood/8, whiting/8, middle/7, low/7, well/7, while/7, tart/7, story/7, song/7, subject/7, chimney/7, lory/7, world/7, pardon/7, puppy/7, oop/7, procession/6, answer/6, lizard/6, youth/6, roof/6, executioner/6, glad/5, watch/5, argument/5, arch/5, case/5, officer/5, morning/5, clock/5, knee/5, teacup/5, back/5, corner/5, sky/5, finger/5, friend/5, hookah/5, flamingo/5, grass/4, fire/4, number/4, age/4, thimble/4, week/4, master/4, fancy/4, sigh/4, 
tag 4:
	,/2419, !/335, '/197, ?/126, and/47, than/11, sob/4, since/4, grave/3, curtsey/2, edition/2, shore/1, bright/1, 
tag 5:
	go/179, <unk>/161, get/108, look/106, come/78, make/76, very/53, turn/42, tell/42, take/39, all/37, not/37, sit/36, as/34, put/34, grow/32, talk/31, give/30, leave/27, too/26, rather/25, so/25, call/25, run/22, quite/21, write/21, repeat/20, keep/20, walk/19, hold/17, find/17, open/17, finish/17, eat/16, nothing/16, fall/15, let/15, live/15, read/14, once/13, shake/13, lie/13, follow/13, swim/12, manage/12, learn/12, matter/11, move/11, play/11, cry/11, interrupt/11, explain/11, please/10, hand/10, help/10, stop/10, sneeze/10, jump/10, drink/10, consider/10, break/10, sing/10, order/9, bring/9, hurry/9, sigh/8, throw/8, watch/8, become/8, drop/7, fetch/7, shut/7, cross/6, miss/6, trouble/6, mark/6, stay/6, deny/6, nurse/6, carry/6, fly/5, choke/5, wave/5, rise/5, dig/5, fast/5, reach/5, shout/5, lead/5, puzzle/5, attend/5, cut/5, set/5, point/5, settle/5, mutter/5, lose/5, bow/5, send/5, show/5, nibble/5, 
tag 6:
	<unk>/262, thing/80, time/74, turtle/60, way/57, voice/51, rabbit/50, one/44, tone/42, day/33, minute/32, hare/31, cat/30, sort/23, child/21, question/21, side/21, end/19, foot/18, remark/17, soup/16, garden/16, pig/16, high/16, use/15, bit/15, idea/15, creature/14, size/14, silence/14, book/13, room/13, bird/12, serpent/12, fan/12, rest/12, deal/12, word/11, conversation/11, ear/11, surprise/11, hurry/11, glove/11, sister/11, place/11, witness/10, hall/9, tea/9, fish/9, party/9, top/9, dream/9, moral/9, key/9, piece/9, nose/8, mushroom/8, sound/8, verse/8, gardener/8, opportunity/8, hair/7, puzzle/7, bat/7, people/7, pepper/7, bread/7, rule/7, temper/7, business/7, ground/7, pocket/7, shriek/7, chin/7, shoe/7, neck/7, paw/7, crowd/7, girl/7, kind/7, sentence/6, egg/6, dog/6, English/6, evidence/6, water/6, politely/6, thought/6, adventure/6, smile/6, cake/6, animal/6, man/6, dish/6, pair/6, fear/5, judge/5, hour/5, juror/5, likely/5, meaning/5, 
tag 7:
	be/411, have/274, n't/217, to/159, do/134, will/99, not/79, could/78, would/75, must/44, seem/40, never/39, all/37, should/32, only/30, can/30, ca/28, might/28, begin/26, quite/21, shall/20, just/20, ever/17, ought/14, may/14, soon/13, wo/12, hardly/12, always/11, even/9, like/7, certainly/7, well/7, venture/7, keep/6, hastily/6, almost/6, remember/6, instantly/5, dare/5, get/5, suddenly/4, ready/4, both/4, simply/3, decide/3, sha/3, immediately/3, really/3, oblige/3, possibly/3, Sha/3, need/3, six/2, dip/2, cautiously/2, nowhere/2, still/2, usually/2, mine/2, thoroughly/2, number/1, cat/1, promise/1, angry/1, next/1, 
tag 8:
	on/127, out/117, it/115, up/100, down/86, off/73, you/59, <unk>/59, me/54, much/51, about/47, over/40, them/39, back/35, her/32, more/31, well/29, so/29, him/29, round/27, away/25, one/23, in/22, herself/22, enough/18, to/17, something/16, close/14, itself/14, mad/14, us/14, half/13, large/13, far/13, anxiously/13, soon/12, with/11, long/11, dry/10, yourself/10, some/10, near/10, else/10, small/10, timidly/9, any/9, from/9, together/9, angrily/9, lesson/9, late/8, slowly/8, asleep/8, suddenly/7, silent/7, croquet/7, myself/7, right/7, taste/6, through/6, nonsense/6, frighten/6, present/6, hard/6, home/5, these/5, sadly/5, nearly/5, across/5, tired/5, different/5, wrong/5, somebody/5, upset/4, true/4, thoughtfully/4, care/4, splash/4, ready/4, glass/4, violently/4, solemnly/4, severely/4, people/4, deeply/4, alone/4, pale/4, many/3, everything/3, trot/3, bend/3, whatever/3, loudly/3, hunt/3, instead/3, em/3, stare/3, dull/3, certain/3, themselves/3, below/3, 
tag 9:
	she/553, i/544, you/321, it/300, alice/169, they/152, he/125, there/70, who/55, to/38, that/35, what/35, we/34, which/26, this/23, then/17, how/16, do/13, wo/12, without/8, would/8, everybody/8, nobody/5, often/5, half/4, seven/4, one/3, nothing/3, nearly/3, mustard/3, soup/2, something/2, latitude/2, else/2, somebody/2, merely/2, hurriedly/2, people/2, ferret/2, attempt/1, crumb/1, 
tag 10:
	and/350, say/333, S/103, s/94, with/24, thought/17, cry/15, think/10, after/7, add/7, half/6, shout/6, continue/6, please/5, father/4, yer/4, scream/4, feel/3, '/3, exclaim/3, rub/3, plead/3, everything/2, mostly/2, rule/1, William/1, frown/1, use/1, busily/1, 
tag 11:
	the/1635, a/627, her/162, his/96, this/70, no/67, your/62, my/58, an/57, its/57, their/52, very/51, some/41, one/35, any/30, that/27, another/21, all/17, every/12, those/10, two/9, what/8, four/8, our/8, which/8, each/8, soo/7, these/7, o/6, both/5, catch/4, beau/4, several/4, gently/3, nine/3, personal/2, always/2, whose/2, course/1, quite/1, 
tag 12:
	little/127, <unk>/108, mock/55, good/39, very/39, great/39, march/35, white/30, large/28, other/26, three/26, same/24, first/22, long/22, curious/21, right/20, more/18, poor/17, next/16, last/14, queer/13, whole/13, old/11, foot/11, beautiful/10, low/10, own/10, tea/9, many/9, few/9, deep/8, most/8, stupid/7, hot/7, bright/7, golden/7, cheshire/7, loud/7, dreadfully/6, guinea/6, two/6, glass/6, sharp/6, melancholy/6, inch/6, nice/6, e/6, shrill/5, simple/5, small/5, young/5, only/5, kid/5, offend/5, strange/5, mile/5, tremble/5, sudden/5, new/5, dear/4, ten/4, flower/4, second/4, capital/4, tiny/4, dead/4, short/4, different/4, lesson/4, confusing/3, dark/3, bark/3, croquet/3, sleepy/3, real/3, really/3, sulky/3, bad/3, grand/3, unfortunate/3, proper/3, hurried/3, timid/3, lovely/3, solemn/3, general/3, red/3, complain/2, clever/2, never/2, wonder/2, open/2, railway/2, quiet/2, coax/2, parchment/2, important/2, feeble/2, crimson/2, encouraging/2, least/2, 
tag 13:
	to/515, of/513, in/345, at/212, and/176, with/146, for/107, into/67, on/66, by/58, about/47, like/43, or/41, after/36, all/34, that/29, from/26, upon/26, such/24, as/20, take/18, without/18, under/16, round/14, than/13, behind/13, quite/12, among/12, poor/10, near/10, against/9, here/8, through/8, eat/6, along/6, between/6, both/5, shrink/5, except/4, blow/4, toss/4, past/3, above/3, raise/3, examine/3, around/3, nine/2, beautifully/2, heavy/2, hide/2, within/2, 
tag 14:
	be/530, say/199, think/107, know/107, see/96, have/70, begin/66, do/64, try/45, hear/45, like/42, find/40, speak/37, reply/34, feel/34, ask/33, wonder/24, mean/23, change/22, happen/21, wish/21, wait/20, such/17, add/17, want/15, stand/14, remember/14, suppose/14, forget/13, afraid/12, join/12, draw/12, appear/11, notice/11, listen/10, believe/10, remark/9, beg/9, could/8, mind/8, execute/8, answer/8, understand/7, grin/7, pass/7, generally/7, box/7, hope/7, offend/6, whisper/6, guess/6, run/6, fancy/6, sound/5, vanish/5, quietly/5, laugh/5, grunt/5, remain/5, prove/5, fit/5, meet/5, catch/5, hole/4, teach/4, worth/4, expect/4, fold/4, growl/3, bite/3, sleep/3, win/3, continue/3, address/3, advance/3, check/3, impossible/3, mention/3, exclaim/3, choose/3, breathe/3, hate/3, tuck/3, cheer/3, bear/3, clear/2, share/2, encourage/2, conclude/2, kindly/2, signify/2, annoy/2, murder/2, advisable/2, pour/2, kneel/2, ornament/2, occur/2, steal/2, 
tag 15:
	and/297, as/209, but/170, that/131, so/97, if/96, what/92, when/79, then/56, how/52, do/47, for/46, or/36, just/32, sure/24, till/21, <unk>/21, while/19, perhaps/17, now/15, because/15, where/15, only/15, would/13, anything/13, before/13, will/12, still/11, though/11, whether/11, either/10, yet/10, all/9, even/9, exactly/8, who/7, can/6, glad/6, why/5, ti/5, which/5, until/5, thank/5, shall/5, label/3, nor/3, suddenly/2, mary/2, stamp/2, luckily/2, wherever/2, unless/2, presently/2, pennyworth/2, ]/2, weak/2, Mary/2, funny/2, longitude/2, nonsense/1, ?/1, 
tag 16:
	again/83, <unk>/74, it/64, that/53, oh/45, now/45, here/43, why/35, well/34, there/29, not/28, before/25, no/23, then/21, however/20, first/18, come/18, indeed/16, dear/16, bill/16, yet/15, down/15, dinah/14, next/13, yes/13, please/11, old/11, hastily/10, twinkle/9, two/8, eagerly/8, five/8, never/7, sir/7, really/7, let/7, certainly/7, prise/6, wow/6, altogether/5, impatiently/5, usual/5, important/5, ah/5, aloud/5, unimportant/5, decidedly/4, indignantly/4, outside/4, sharply/4, behead/4, stuff/4, pat/4, yawn/4, besides/4, alas/4, otherwise/4, beautiful/3, ma/3, frown/3, music/3, pant/3, tortoise/3, riddle/3, stair/3, growl/3, love/3, quadrille/3, pray/3, fury/3, hush/3, Ann/3, flamingo/2, twelve/2, treacle/2, uglification/2, remove/2, three/2, doubtfully/2, already/2, chain/2, mystery/2, goose/2, sh/2, rude/2, sooner/2, triumphantly/2, ugh/2, patiently/2, thump/2, morcar/2, tut/2, inside/2, French/2, contemptuously/2, Edwin/2, useful/2, afterwards/1, 
tag 17:
	alice/225, it/116, herself/61, her/54, them/49, all/47, <unk>/45, this/41, that/40, you/31, course/25, once/21, last/21, two/17, nothing/15, him/14, me/14, sight/10, tear/10, which/10, heart/10, everything/9, work/9, do/9, first/9, mine/8, least/7, butter/6, what/6, himself/6, school/6, fact/6, right/5, delight/5, sob/5, knock/5, twice/5, William/5, execution/5, ever/5, sometimes/5, anything/4, particular/4, green/4, honour/4, confusion/4, none/4, hers/4, nobody/3, curl/3, yours/3, apple/3, hearing/3, hungry/3, passion/3, mabel/3, livery/3, whisker/3, instance/3, card/3, mark/2, angry/2, giddy/2, directly/2, shilling/2, flat/2, knot/2, whistle/2, couple/2, wag/2, sell/2, law/2, uglify/2, rapidly/2, tiptoe/2, wonderland/2, secondly/2, custody/2, furrow/2, paris/2, front/2, mercia/2, william/2, argue/2, northumbria/2, fright/2, pretend/2, dinn/2, another/1, not/1, cupboard/1, 
tag 18:
	use/17, set/10, rattle/3, alice/2, care/2, hard/2, lap/2, picture/1, alarm/1, sneeze/1, share/1, 
tag 19:
	these/2, ten/2, first/2, afterwards/1, verdict/1, hour/1, <unk>/1, carry/1, come/1, soldier/1, sentence/1, 
```
