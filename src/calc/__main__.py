"""Read-only command-line discovery for the algebra/calculus atlas."""

import argparse
import json
from .catalog import TOPICS, get_topic, list_topics
from .evidence import get_evidence


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="calc-atlas", description="Explore calculus and algebra topics")
    sub = parser.add_subparsers(dest="command", required=True)
    listing = sub.add_parser("topics", help="list indexed topics")
    listing.add_argument("--family", default=None)
    detail = sub.add_parser("topic", help="show one topic")
    detail.add_argument("slug")
    evidence = sub.add_parser("evidence", help="show source and test evidence for a topic")
    evidence.add_argument("slug")
    evidence.add_argument("--json", action="store_true", help="machine-readable record")
    args = parser.parse_args(argv)
    if args.command == "topics":
        topics = list_topics(args.family)
        if args.family is not None and not topics:
            parser.error(f"unknown family: {args.family}")
        for topic in topics:
            print(f"{topic.slug:36} {topic.status:11} {topic.family}")
        print(f"{len(topics)} topics")
        return 0
    if args.command == "evidence":
        try:
            record = get_evidence(args.slug)
        except KeyError:
            parser.error(f"unknown topic: {args.slug}")
        if args.json:
            print(json.dumps(record, indent=2, ensure_ascii=False, sort_keys=True))
        else:
            print(f"{record['slug']}\nstatus: {record['status']}"
                  f"\nevidence: {record['evidence_level']}"
                  f"\nscope: {record['scope']}"
                  f"\nguide: {record['guide'] or 'not written'}"
                  f"\nsources: {', '.join(record['source_ids'])}"
                  f"\nformal proof: {record['formal_proof']}")
        return 0
    try:
        topic = get_topic(args.slug)
    except KeyError:
        parser.error(f"unknown topic: {args.slug}")
    print(f"{topic.title}\nslug: {topic.slug}\nfamily: {topic.family}\nstatus: {topic.status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

