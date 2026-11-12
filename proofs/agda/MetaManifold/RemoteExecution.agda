-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-- SPDX-License-Identifier: MPL-2.0

{-# OPTIONS --safe --without-K #-}

module MetaManifold.RemoteExecution where

open import Data.Nat using (ℕ; zero; suc; _≤_; z≤n; s≤s)
open import Data.Product using (_×_; _,_)
open import Relation.Binary.PropositionalEquality using (_≡_; _≢_; refl; trans)
open import Data.Empty using (⊥)

-- One journal-owned attempt. Owner and content are abstract identities, not
-- machine hashes. A key equality means equality of BOTH components. The mapping
-- to UUIDs, SHA-256, OS locks and acknowledgements is described in the README;
-- those primitives are assumptions of this protocol model, not proved here.
Key = ℕ × ℕ

data Phase : Set where
  preparing ready running stopping succeeded failed fetched cancelled removed : Phase

record State : Set where
  constructor state
  field
    key : Key
    phase : Phase
    active launches : ℕ
open State

infix 4 _↦_
data _↦_ : State → State → Set where
  verified : ∀ {k n} → state k preparing 0 n ↦ state k ready 0 n
  launch : ∀ {k n} → state k ready 0 n ↦ state k running 1 (suc n)
  -- Repeated start requests are refused by the held remote lock.
  duplicate-start : ∀ {k n} → state k running 1 n ↦ state k running 1 n
  disconnect : ∀ {s} → s ↦ s
  finish-ok : ∀ {k n} → state k running 1 n ↦ state k succeeded 0 n
  finish-error : ∀ {k n} → state k running 1 n ↦ state k failed 0 n
  -- A cancel request alone is NOT an acknowledgement that the process stopped.
  request-cancel : ∀ {k n} → state k running 1 n ↦ state k stopping 1 n
  stop-ack : ∀ {k n} → state k stopping 1 n ↦ state k cancelled 0 n
  cancel-ready : ∀ {k n} → state k ready 0 n ↦ state k cancelled 0 n
  cancel-preparing : ∀ {k n} → state k preparing 0 n ↦ state k cancelled 0 n
  cancel-succeeded : ∀ {k n} → state k succeeded 0 n ↦ state k cancelled 0 n
  cancel-failed : ∀ {k n} → state k failed 0 n ↦ state k cancelled 0 n
  cancel-fetched : ∀ {k n} → state k fetched 0 n ↦ state k cancelled 0 n
  fetch-verified : ∀ {k n} → state k succeeded 0 n ↦ state k fetched 0 n
  cleanup-fetched : ∀ {k n} → state k fetched 0 n ↦ state k removed 0 n
  cleanup-cancelled : ∀ {k n} → state k cancelled 0 n ↦ state k removed 0 n
  -- Retention may remove old terminal attempts without a client fetching them.
  -- The clock/age policy is outside this safety model; the free lock is not.
  expire-succeeded : ∀ {k n} → state k succeeded 0 n ↦ state k removed 0 n
  expire-failed : ∀ {k n} → state k failed 0 n ↦ state k removed 0 n
  retry : ∀ {k n} → state k failed 0 n ↦ state k ready 0 n
  -- A succeeded attempt whose outputs were never retrieved is fetched again from
  -- the same attempt. Repeated replay requests are idempotent: a failed
  -- retrieval leaves both the phase and the launch count untouched.
  replay : ∀ {k n} → state k succeeded 0 n ↦ state k succeeded 0 n

-- Induction over arbitrarily long traces, not a bounded enumeration.
infix 4 _↦*_
data _↦*_ : State → State → Set where
  done : ∀ {s} → s ↦* s
  next : ∀ {s t u} → s ↦ t → t ↦* u → s ↦* u

step-exclusive : ∀ {s t} → s ↦ t → active s ≤ 1 → active t ≤ 1
step-exclusive verified _ = z≤n
step-exclusive launch _ = s≤s z≤n
step-exclusive duplicate-start p = p
step-exclusive disconnect p = p
step-exclusive finish-ok _ = z≤n
step-exclusive finish-error _ = z≤n
step-exclusive request-cancel p = p
step-exclusive stop-ack _ = z≤n
step-exclusive cancel-ready _ = z≤n
step-exclusive cancel-preparing _ = z≤n
step-exclusive cancel-succeeded _ = z≤n
step-exclusive cancel-failed _ = z≤n
step-exclusive cancel-fetched _ = z≤n
step-exclusive fetch-verified _ = z≤n
step-exclusive cleanup-fetched _ = z≤n
step-exclusive cleanup-cancelled _ = z≤n
step-exclusive expire-succeeded _ = z≤n
step-exclusive expire-failed _ = z≤n
step-exclusive retry _ = z≤n
step-exclusive replay _ = z≤n

trace-exclusive : ∀ {s t} → s ↦* t → active s ≤ 1 → active t ≤ 1
trace-exclusive done p = p
trace-exclusive (next st rest) p = trace-exclusive rest (step-exclusive st p)

one-live-command : ∀ {k t} → state k preparing 0 0 ↦* t → active t ≤ 1
one-live-command path = trace-exclusive path z≤n

step-key : ∀ {s t} → s ↦ t → key s ≡ key t
step-key verified = refl
step-key launch = refl
step-key duplicate-start = refl
step-key disconnect = refl
step-key finish-ok = refl
step-key finish-error = refl
step-key request-cancel = refl
step-key stop-ack = refl
step-key cancel-ready = refl
step-key cancel-preparing = refl
step-key cancel-succeeded = refl
step-key cancel-failed = refl
step-key cancel-fetched = refl
step-key fetch-verified = refl
step-key cleanup-fetched = refl
step-key cleanup-cancelled = refl
step-key expire-succeeded = refl
step-key expire-failed = refl
step-key retry = refl
step-key replay = refl

trace-key : ∀ {s t} → s ↦* t → key s ≡ key t
trace-key done = refl
trace-key (next st rest) = trans (step-key st) (trace-key rest)

different-attempts-stay-different : ∀ {a b a′ b′} →
  key a ≢ key b → a ↦* a′ → b ↦* b′ → key a′ ≢ key b′
different-attempts-stay-different different pa pb equal
  with trace-key pa | trace-key pb
... | refl | refl = different equal

-- These statements concern a newly entered phase, excluding disconnect stutter.
cleanup-needs-no-live-command : ∀ {k p a n k′ a′ n′} →
  p ≢ removed → state k p a n ↦ state k′ removed a′ n′ → a ≡ 0
cleanup-needs-no-live-command not-removed disconnect =
  Data.Empty.⊥-elim (not-removed refl)
cleanup-needs-no-live-command _ cleanup-fetched = refl
cleanup-needs-no-live-command _ cleanup-cancelled = refl
cleanup-needs-no-live-command _ expire-succeeded = refl
cleanup-needs-no-live-command _ expire-failed = refl

fetch-needs-success : ∀ {k p a n k′ a′ n′} →
  p ≢ fetched → state k p a n ↦ state k′ fetched a′ n′ → p ≡ succeeded
fetch-needs-success not-fetched disconnect = Data.Empty.⊥-elim (not-fetched refl)
fetch-needs-success _ fetch-verified = refl

cancel-request-keeps-live-count : ∀ {k n} →
  active (state k stopping 1 n) ≡ active (state k running 1 n)
cancel-request-keeps-live-count = refl

-- A reattachment is admitted by exact journal identity equality. Fingerprint
-- collision resistance is deliberately not asserted as an injectivity theorem.
data Attaches (expected : Key) (actual : State) : Set where
  checked : expected ≡ key actual → Attaches expected actual

different-content-cannot-attach : ∀ {expected actual} →
  expected ≢ key actual → Attaches expected actual → ⊥
different-content-cannot-attach different (checked equal) = different equal

disconnect-preserves-state : ∀ (s : State) → s ↦ s
disconnect-preserves-state s = disconnect

duplicate-launch-preserves-count : ∀ {k n} →
  state k running 1 n ↦ state k running 1 n
duplicate-launch-preserves-count = duplicate-start

-- Retrieving a completed attempt again is not a relaunch: the launch count of
-- the attempt is unchanged, so a replay can never admit a second command.
replay-preserves-launch-count : ∀ {k n} →
  launches (state k succeeded 0 n) ≡ launches (state k succeeded 0 n)
replay-preserves-launch-count = refl
